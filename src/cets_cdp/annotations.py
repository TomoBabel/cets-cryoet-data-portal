"""Portal annotations ⇄ CETS annotations.

Portal → CETS: every ``AnnotationFile`` of a run binds to the region's tomogram at the same
``tomogram_voxel_spacing_id`` (the portal tomogram the file was made on). Point files (ndjson, voxel indices at
that voxel spacing, corner origin, no half-voxel term — ``common/point_converter.py``) become ``PointSet3D`` /
``PointMatrixSet3D`` in the tomogram's centred physical frame: ``p = (loc − floor(N/2)) · s``;
``xyz_rotation_matrix`` is RELION's ``A(rot, tilt, psi)`` = the particle→tomogram matrix and is used as is.
Masks (OME-Zarr + MRC, ``importers/annotation.py``) become ``SegmentationMask3D`` with the grid read from the
zarr metadata / MRC header, the file referenced by URL.

CETS → portal: point annotations are staged as ``relion4_star`` files (``rlnCoordinateX/Y/Z`` in voxels of the
bound tomogram, ``rlnOrigin*Angst = 0``, Eulers from the matrices, ``rlnImagePixelSize`` = the voxel size) — the
only accepted point source format that carries orientations without a copick dependency; the backend divides by
``binning`` (1) and stores ``inv(scipy ZYZ)`` = the same matrix (``point_converter.py:315-363``). Local masks are
symlinked; remote ones are left to the curator with the URL.
"""

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import starfile
from cryoet_alignment.io.cets.annotations import (
    POINT_KINDS,
    ResolvedPoints,
    annotation_id,
    annotation_kind,
    annotation_points,
    mask_entity,
    point_set_entity,
    tomogram_frame,
)
from cryoet_alignment.io.cets.cli_support import Gate, SeriesReport
from cryoet_alignment.io.cets.companion import AnnotationCompanion
from cryoet_alignment.io.cets.euler import matrices_to_zyz, zyz_to_matrices
from cryoet_alignment.io.cets.rotation import check_rotation

from cets_cdp.api import PortalAnnotation, PortalAnnotationFile, PortalRunData, https_to_s3

POINT_SHAPES = ("Point", "OrientedPoint", "InstanceSegmentation")
MASK_SHAPES = ("SegmentationMask", "InstanceSegmentationMask", "SemanticSegmentationMask")
DEFAULT_SHAPES = POINT_SHAPES + MASK_SHAPES


# --------------------------------------------------------------------------- portal → CETS


@dataclass
class NdjsonPoints:
    locations: np.ndarray  # (N,3) voxel indices
    matrices: Optional[np.ndarray]  # (N,3,3)
    instance_ids: Optional[List[int]]
    kinds: List[str]


def read_ndjson_points(path: Path) -> NdjsonPoints:
    locs: List[List[float]] = []
    mats: List[Any] = []
    inst: List[Any] = []
    kinds: List[str] = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        loc = rec["location"]
        locs.append([float(loc["x"]), float(loc["y"]), float(loc["z"])])
        kinds.append(str(rec.get("type", "point")))
        mats.append(rec.get("xyz_rotation_matrix"))
        inst.append(rec.get("instance_id"))
    if not locs:
        raise ValueError(f"{path}: no point records")
    matrices = None
    if all(m is not None for m in mats):
        matrices = np.asarray(mats, dtype=np.float64).reshape(-1, 3, 3)
    elif any(m is not None for m in mats):
        raise ValueError(f"{path}: some records carry xyz_rotation_matrix and some do not")
    instance_ids = [int(v) for v in inst] if all(v is not None for v in inst) else None
    return NdjsonPoints(np.asarray(locs, dtype=np.float64), matrices, instance_ids, kinds)


@dataclass
class AnnotationImport:
    entities: List[Any] = field(default_factory=list)
    companions: Dict[str, AnnotationCompanion] = field(default_factory=dict)


def _bind_tomogram(data: PortalRunData, f: PortalAnnotationFile, tomograms_by_portal_id: Dict[int, Any]):
    """The region's tomogram entity the file is attached to (same portal voxel-spacing id; portal standard first)."""
    cands = [
        tid
        for tid, tvs in data.tomogram_voxel_spacing_ids.items()
        if tvs == f.tomogram_voxel_spacing_id and tid in tomograms_by_portal_id
    ]
    if not cands:
        # older recordings without voxel-spacing ids: match by voxel spacing
        cands = [
            t.id
            for t in data.tomograms
            if f.voxel_spacing is not None
            and abs(t.voxel_spacing - f.voxel_spacing) < 1e-3
            and t.id in tomograms_by_portal_id
        ]
    if not cands:
        return None
    ranked = sorted(
        cands, key=lambda tid: (not bool(next(t.is_portal_standard for t in data.tomograms if t.id == tid)), tid)
    )
    return tomograms_by_portal_id[ranked[0]]


def _companion_base(a: PortalAnnotation, f: PortalAnnotationFile, tomo_id: str, kind: str) -> AnnotationCompanion:
    return AnnotationCompanion(
        kind=kind,
        tomogram_id=tomo_id,
        source_tool="cryoET Data Portal",
        source_ref=f"portal annotation {a.id} file {f.id} ({f.shape_type}, {f.format})",
        name=a.object_name,
        portal_annotation_id=a.id,
        object_name=a.object_name,
        object_id=a.object_id,
        object_state=a.object_state,
        annotation_method=a.annotation_method,
        method_type=a.method_type,
        ground_truth_status=a.ground_truth_status,
        is_curator_recommended=a.is_curator_recommended,
        annotation_software=a.annotation_software,
        confidence_precision=a.confidence_precision,
        confidence_recall=a.confidence_recall,
        metadata=a.metadata,
        shape_type=f.shape_type,
        file_format=f.format,
        file_path=f.https_path,
        voxel_spacing_a=f.voxel_spacing,
        is_visualization_default=f.is_visualization_default,
    )


def portal_annotations_to_cets(
    data: PortalRunData,
    tomograms_by_portal_id: Dict[int, Any],
    sr: SeriesReport,
    *,
    alignment_id: Optional[int],
    shapes: Sequence[str] = DEFAULT_SHAPES,
    scheme: str = "https",
) -> AnnotationImport:
    """Convert the run's annotation files to CETS entities bound to the region's tomograms (keyed by portal id)."""
    out = AnnotationImport()
    by_portal_id = tomograms_by_portal_id
    n_points_ok = 0
    n_points = 0
    n_grid_ok = 0
    n_masks = 0
    for a in data.annotations:
        for f in a.files:
            if f.shape_type not in shapes:
                continue
            if f.shape_type in MASK_SHAPES and f.format != "zarr":
                continue  # the MRC twin is recorded on the zarr's companion entry
            if f.shape_type in POINT_SHAPES and f.format != "ndjson":
                continue
            if alignment_id is not None and f.alignment_id is not None and f.alignment_id != alignment_id:
                sr.warnings.append(
                    f"annotation {a.id} ({a.object_name}, {f.shape_type}) belongs to alignment {f.alignment_id}, "
                    f"not the exported {alignment_id}: skipped",
                )
                continue
            tomo = _bind_tomogram(data, f, by_portal_id)
            if tomo is None:
                sr.warnings.append(
                    f"annotation {a.id} ({a.object_name}, {f.shape_type}) is attached to voxel spacing "
                    f"{f.voxel_spacing} Å which has no tomogram in the region: skipped",
                )
                continue
            frame = tomogram_frame(tomo)
            key = f"{a.id}_{f.shape_type.lower()}"
            ann_id = annotation_id(str(tomo.id), key)
            if f.shape_type in POINT_SHAPES:
                if not f.local_path:
                    sr.warnings.append(f"annotation {a.id} ({f.shape_type}): ndjson not cached: skipped")
                    continue
                pts = read_ndjson_points(Path(f.local_path))
                if f.voxel_spacing is not None and abs(f.voxel_spacing - frame.spacing_a) > 1e-3:
                    sr.warnings.append(
                        f"annotation {a.id}: file voxel spacing {f.voxel_spacing} differs from tomogram {tomo.id} ({frame.spacing_a})",
                    )
                corner = pts.locations * frame.spacing_a
                inside = frame.inside(corner)
                n_points += len(corner)
                n_points_ok += int(inside.sum())
                mats = None
                if pts.matrices is not None:
                    for i, r in enumerate(pts.matrices):
                        check_rotation(r, f"annotation {a.id} record {i} xyz_rotation_matrix")
                    mats = pts.matrices
                ent = point_set_entity(
                    annotation_id=ann_id,
                    tomogram_id=str(tomo.id),
                    points_a=frame.corner_to_cets(corner),
                    matrices=mats,
                    name=a.object_name,
                )
                comp = _companion_base(a, f, str(tomo.id), "oriented_points" if mats is not None else "points")
                if pts.instance_ids is not None:
                    comp.instance_ids = pts.instance_ids
                    comp.kind = "instance_points"
                out.entities.append(ent)
                out.companions[ann_id] = comp
            else:
                grid = f.grid
                voxel = f.grid_voxel_a or f.voxel_spacing
                if grid is None:
                    sr.warnings.append(
                        f"annotation {a.id} ({a.object_name}, mask): grid not readable from the zarr/mrc metadata; "
                        f"assuming the tomogram grid {frame.size_px}",
                    )
                    grid = frame.size_px
                if voxel is None:
                    voxel = frame.spacing_a
                n_masks += 1
                if tuple(grid) == tuple(frame.size_px) and abs(float(voxel) - frame.spacing_a) < 1e-3:
                    n_grid_ok += 1
                else:
                    sr.warnings.append(
                        f"annotation {a.id} ({a.object_name}, mask): grid {tuple(grid)} @ {voxel} Å differs from tomogram "
                        f"{tomo.id} {frame.size_px} @ {frame.spacing_a} Å; the mask keeps its own frame",
                    )
                ent = mask_entity(
                    annotation_id=ann_id,
                    tomogram_id=str(tomo.id),
                    path=https_to_s3(f.https_path) if scheme == "s3" else f.https_path,
                    size_px=grid,
                    voxel_size_a=float(voxel),
                    name=a.object_name,
                )
                comp = _companion_base(a, f, str(tomo.id), "mask")
                mrc = next((x for x in a.files if x.shape_type == f.shape_type and x.format == "mrc"), None)
                if mrc is not None:
                    comp.mrc_path = https_to_s3(mrc.https_path) if scheme == "s3" else mrc.https_path
                meta = a.metadata or {}
                comp.mask_label = 1  # portal masks are binarised at ingestion (importers/annotation.py:330-351)
                if f.shape_type == "InstanceSegmentationMask":
                    comp.mask_label = None
                    comp.kind = "instance_mask"
                comp.name = a.object_name or meta.get("annotation_object", {}).get("name")
                out.entities.append(ent)
                out.companions[ann_id] = comp
    if n_points and n_points_ok != n_points:
        sr.warnings.append(f"{n_points - n_points_ok} of {n_points} annotation points lie outside their tomogram grid")
    if n_masks:
        sr.gates.append(Gate("mask_grid_matches_tomogram", n_grid_ok == n_masks, value=n_grid_ok, expected=n_masks))
    return out


# --------------------------------------------------------------------------- CETS → portal staging


@dataclass
class StagedAnnotation:
    annotation_id: str
    kind: str  # points | oriented_points | mask
    object_name: Optional[str]
    metadata: Optional[dict]
    glob_rel: Optional[str]  # staging-relative path of the staged file (None when nothing could be staged)
    shape_key: str  # Point | OrientedPoint | SemanticSegmentationMask
    file_format: str  # relion4_star | zarr | mrc
    is_visualization_default: Optional[bool]
    mask_label: Optional[int]
    note: Optional[str] = None
    n: int = 0


def stage_points_star(
    resolved: ResolvedPoints, path: Path, *, run: str, extra_columns: Optional[Dict[str, list]] = None
) -> Path:
    """Write a ``relion4_star`` file the backend's ``from_relion4_star`` turns back into the same voxel locations
    (``rlnTomoName`` = the run name so a ``filter_value: '{run_name}'`` source also works)."""
    frame = resolved.frame
    vox = resolved.points_corner_a / frame.spacing_a
    n = resolved.n
    data: Dict[str, list] = {
        "rlnTomoName": [run] * n,
        "rlnCoordinateX": vox[:, 0].tolist(),
        "rlnCoordinateY": vox[:, 1].tolist(),
        "rlnCoordinateZ": vox[:, 2].tolist(),
        "rlnOriginXAngst": [0.0] * n,
        "rlnOriginYAngst": [0.0] * n,
        "rlnOriginZAngst": [0.0] * n,
        "rlnOpticsGroup": [1] * n,
    }
    if resolved.matrices is not None:
        e = matrices_to_zyz(resolved.matrices)
        data["rlnAngleRot"], data["rlnAngleTilt"], data["rlnAnglePsi"] = (
            e[:, 0].tolist(),
            e[:, 1].tolist(),
            e[:, 2].tolist(),
        )
    else:
        data["rlnAngleRot"] = data["rlnAngleTilt"] = data["rlnAnglePsi"] = [0.0] * n
    for k, v in (extra_columns or {}).items():
        if len(v) == n and k not in data:
            data[k] = list(v)
    optics = pd.DataFrame(
        {"rlnOpticsGroup": [1], "rlnOpticsGroupName": ["opticsGroup1"], "rlnImagePixelSize": [frame.spacing_a]}
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    starfile.write({"optics": optics, "particles": pd.DataFrame(data)}, str(path), float_format="%.6f")
    return path


def relion4_star_to_points(path: Path) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """The backend's ``_from_relion4_star`` arithmetic (point_converter.py:315-363) with binning 1: locations =
    coordinate − originAngst / rlnImagePixelSize; matrices = inv(scipy ZYZ) = ``zyz_to_matrix``."""
    blocks = starfile.read(str(path), always_dict=True)
    p = blocks["particles"]
    pix = float(blocks["optics"]["rlnImagePixelSize"].iloc[0])
    loc = (
        p[["rlnCoordinateX", "rlnCoordinateY", "rlnCoordinateZ"]].to_numpy(float)
        - p[["rlnOriginXAngst", "rlnOriginYAngst", "rlnOriginZAngst"]].to_numpy(float) / pix
    )
    mats = None
    if all(c in p.columns for c in ("rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi")):
        mats = zyz_to_matrices(p[["rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"]].to_numpy(float))
    return loc, mats


def stage_annotations(
    region,
    companion_annotations: Dict[str, AnnotationCompanion],
    sr: SeriesReport,
    *,
    staging: Path,
    doc_dir: Path,
    run: str,
    annotation_ids: Sequence[str] = (),
    overwrite: bool = False,
) -> List[StagedAnnotation]:
    """Stage every point annotation as a relion4 star and every local mask as a symlink; report the rest."""
    from cets_cdp.from_cets import _link, _resolve_doc_path

    out: List[StagedAnnotation] = []
    anns = list(region.annotations or [])
    if annotation_ids:
        wanted = set(annotation_ids)
        anns = [a for a in anns if str(a.id) in wanted]
    worst = 0.0
    worst_rot = 0.0
    for ann in anns:
        kind = annotation_kind(ann)
        comp = companion_annotations.get(str(ann.id))
        meta = comp.metadata if comp else None
        name = (comp.object_name if comp else None) or getattr(ann, "name", None)
        cmeta = config_metadata(meta)
        if kind in POINT_KINDS:
            r = annotation_points(ann, region)
            shape_key = "OrientedPoint" if r.matrices is not None else "Point"
            d = staging / "annotations" / run
            path = d / f"{file_key(meta, str(ann.id))}_{shape_key.lower()}.star"
            if path.exists() and not overwrite:
                raise FileExistsError(f"{path} exists (use --overwrite)")
            extra = {"rlnTomoParticleId": comp.instance_ids} if comp and comp.instance_ids else None
            stage_points_star(r, path, run=run, extra_columns=extra)
            loc, mats = relion4_star_to_points(path)
            worst = max(worst, float(np.abs(loc - r.points_corner_a / r.frame.spacing_a).max()))
            if r.matrices is not None and mats is not None:
                rel = np.einsum("nij,nkj->nik", mats, r.matrices)
                worst_rot = max(worst_rot, float(np.abs(rel - np.eye(3)).max()))
            out.append(
                StagedAnnotation(
                    annotation_id=str(ann.id),
                    kind=kind,
                    object_name=name,
                    metadata=cmeta,
                    glob_rel=f"annotations/{{run_name}}/{path.name}",
                    shape_key=shape_key,
                    file_format="relion4_star",
                    is_visualization_default=comp.is_visualization_default if comp else None,
                    mask_label=None,
                    n=r.n,
                )
            )
        elif kind == "mask":
            src = _resolve_doc_path(getattr(ann, "path", None), doc_dir)
            remote = str(getattr(ann, "path", "") or "").startswith(("http", "s3:"))
            if src is not None and not remote and src.exists():
                d = staging / "annotations" / run
                d.mkdir(parents=True, exist_ok=True)
                dst = d / f"{file_key(meta, str(ann.id))}_segmentationmask{src.suffix}"
                _link(src, dst, overwrite)
                fmt = "zarr" if src.suffix == ".zarr" or src.is_dir() else "mrc"
                out.append(
                    StagedAnnotation(
                        annotation_id=str(ann.id),
                        kind=kind,
                        object_name=name,
                        metadata=cmeta,
                        glob_rel=f"annotations/{{run_name}}/{dst.name}",
                        shape_key="SemanticSegmentationMask",
                        file_format=fmt,
                        is_visualization_default=comp.is_visualization_default if comp else None,
                        mask_label=comp.mask_label if comp and comp.mask_label is not None else 1,
                    )
                )
            else:
                out.append(
                    StagedAnnotation(
                        annotation_id=str(ann.id),
                        kind=kind,
                        object_name=name,
                        metadata=cmeta,
                        glob_rel=None,
                        shape_key="SemanticSegmentationMask",
                        file_format="zarr",
                        is_visualization_default=comp.is_visualization_default if comp else None,
                        mask_label=comp.mask_label if comp and comp.mask_label is not None else 1,
                        note=f"mask file is not local ({getattr(ann, 'path', None)}): source left as TODO",
                    )
                )
                sr.warnings.append(f"annotation {ann.id}: {out[-1].note}")
        else:
            sr.warnings.append(f"annotation {ann.id} ({kind}) has no portal staging: skipped")
    if any(s.kind in POINT_KINDS for s in out):
        sr.gates.append(
            Gate(
                "backend_point_parser_reproduces_points",
                worst < 1e-6,
                value=worst,
                expected="< 1e-6 voxel (relion4_star arithmetic, binning 1)",
            )
        )
        if worst_rot > 0:
            sr.gates.append(
                Gate(
                    "backend_point_parser_reproduces_rotations",
                    worst_rot < 1e-6,
                    value=worst_rot,
                    expected="< 1e-6 (|R_out R_in^T - I|)",
                )
            )
    return out


#: Keys of the portal's annotation metadata json that are ingestion OUTPUTS, not config inputs.
_METADATA_OUTPUT_KEYS = ("files", "last_updated_at", "deposition_id", "alignment_metadata_path", "object_count")


def config_metadata(meta: Optional[dict]) -> Optional[dict]:
    """The portal metadata json reduced to the ingestion config's ``annotations[].metadata`` fields."""
    if not meta:
        return None
    return {k: v for k, v in meta.items() if k not in _METADATA_OUTPUT_KEYS}


def file_key(meta: Optional[dict], annotation_id: str) -> str:
    """Staged file stem shared by every run of a deposition: ``<deposition>_<ingest id>`` when known (the ingest
    id is unique only within a deposition's config), else the annotation id."""
    ingest = (meta or {}).get("annotation_ingest_id")
    if not ingest:
        return str(annotation_id)
    dep = (meta or {}).get("deposition_id")
    stem = f"{dep}_{ingest}" if dep else str(ingest)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", stem)


PLACEHOLDER_METADATA = {
    "annotation_ingest_id": "TODO(curator)",
    "annotation_object": {"id": "TODO(curator)", "name": "TODO(curator)"},
    "annotation_method": "TODO(curator)",
    "method_type": "TODO(curator)",
    "authors": "TODO(curator)",
    "dates": "TODO(curator)",
    "version": 1.0,
}


def annotation_blocks(staged_by_run: Dict[str, List[StagedAnnotation]], template: Optional[dict]) -> List[dict]:
    """Group staged annotations across runs into ingestion ``annotations`` blocks: annotations with identical
    metadata (the ingestion block is per deposition, identical for every run) share one block per shape key;
    without metadata, by object name."""
    groups: Dict[Tuple[str, str], List[Tuple[str, StagedAnnotation]]] = {}
    for run, items in staged_by_run.items():
        for s in items:
            key = json.dumps(s.metadata, sort_keys=True) if s.metadata else (s.object_name or s.annotation_id)
            groups.setdefault((key, s.shape_key), []).append((run, s))
    tpl_blocks = (template or {}).get("annotations") or []

    def _template_for(obj: Optional[str]) -> Optional[dict]:
        for b in tpl_blocks:
            name = ((b.get("metadata") or {}).get("annotation_object") or {}).get("name")
            if obj and name and name.lower() == obj.lower():
                return b.get("metadata")
        return None

    blocks: List[dict] = []
    for (_key, shape), group in groups.items():
        first = group[0][1]
        obj = first.object_name
        meta = first.metadata or _template_for(obj)
        if meta is None:
            meta = json.loads(json.dumps(PLACEHOLDER_METADATA))
            if obj:
                meta["annotation_object"]["name"] = obj
        globs = sorted({sa.glob_rel for _, sa in group if sa.glob_rel})
        source: Dict[str, Any] = {}
        if shape == "Point":  # Point sources take `columns`, OrientedPoint sources take `order` (schema metadata.yaml)
            source = {"file_format": "relion4_star", "binning": 1, "columns": "xyz"}
        elif shape == "OrientedPoint":
            source = {"file_format": "relion4_star", "binning": 1, "order": "xyz"}
        else:
            source = {
                "file_format": first.file_format,
                "mask_label": first.mask_label if first.mask_label is not None else 1,
            }
        if len(globs) == 1:
            source["glob_string"] = globs[0]
        elif globs:
            source["glob_strings"] = globs
        else:
            source["glob_string"] = "TODO(curator)"
        if first.is_visualization_default is not None:
            source["is_visualization_default"] = bool(first.is_visualization_default)
        blocks.append({"metadata": meta, "sources": [{shape: source}]})
    return blocks
