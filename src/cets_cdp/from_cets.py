"""CETS dataset -> portal staging directory + ingestion-config draft + three separate validations.

Staging (the shape current ingestion configs glob)::

    staging/
      alignment/<run>/<run>.aln                 (all volume_x_rotation == 0 -> format ARETOMO3)
      alignment/<run>/<run>.{xf,tlt,xtilt}      (otherwise -> format IMOD)
      rawtlt/<run>/<run>.rawtlt                 nominal_tilt_angle per raw section
      ctf/<run>/<run>_CTF.txt                   when every image carries CTF metadata (CTFFIND-parsable)
      tiltseries/<run>/<run>.mrc                symlink when the CETS path is local
      tomograms/<run>/<run>.<mrc|zarr>          symlink when the CETS tomogram path is local
      run_to_data_map.tsv                       per-run values when they differ between runs
      ingestion_config.yaml                     the draft
      cets_cdp.report.json

Nothing is inferred: ``processing``, ``reconstruction_method`` etc. come from the companion or the
``--template`` config, else stay ``TODO(curator)`` placeholders.
"""

import copy
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import yaml
from cryoet_alignment.io.aretomo3 import AreTomo3ALN, AreTomo3CTF
from cryoet_alignment.io.cets import ctf as cets_ctf
from cryoet_alignment.io.cets.alignment import (
    ReferenceVolume,
    alignment_from_cets,
    alignment_name_of,
    select_alignment,
    select_tomogram,
)
from cryoet_alignment.io.cets.cli_support import Gate, SeriesReport
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.config import Resolver
from cryoet_alignment.io.cets.frames import FRAME_CONVENTIONS, find_by_id, image_frame
from cryoet_alignment.io.cryoet_data_portal import Alignment
from cryoet_alignment.io.imod import ImodAlignment, ImodTLT, ImodXF, ImodXTILT

TODO = "TODO(curator)"
CONFIG_VERSION = "1.1.0"
METHOD_TYPES = ("fiducial_based", "patch_tracking", "projection_matching", "simulated", "undefined")


# ------------------------------------------------------------------------------------------ per run


def stage_run(
    region,
    res: Resolver,
    sr: SeriesReport,
    *,
    staging: Path,
    doc_dir: Path,
    companion: Optional[Companion],
    alignment_selector=None,
    tomogram_selector: Optional[str] = None,
    overwrite: bool = False,
) -> Dict[str, Any]:
    """Stage one region; returns the per-run values the config needs."""
    cets_alignment = select_alignment(region, alignment_selector)
    ts = find_by_id(region.tilt_series, cets_alignment.tilt_series_id, "tilt series")
    run = ts.id
    aln_name = alignment_name_of(cets_alignment)
    comp_aln = companion.alignment(ts.id, aln_name) if companion else None
    comp_ts = companion.tilt_series.get(ts.id) if companion else None
    comp_images = comp_ts.images if comp_ts else {}

    tomo = select_tomogram(region, cets_alignment, tomogram_selector, companion)
    reference = ReferenceVolume.from_tomogram(tomo)
    native = comp_aln.native_volume_dimension_a if comp_aln and comp_aln.native_volume_dimension_a else None
    hub: Alignment = alignment_from_cets(
        cets_alignment,
        tilt_series=ts,
        reference=reference,
        target_frame=FRAME_CONVENTIONS["ARETOMO3"],
        native_dimension_a=native,
        format_="ARETOMO3",
    )
    res.resolve(
        "reference_tomogram",
        discovered=tomo.id,
        note="companion" if comp_aln and comp_aln.reference_tomogram_id else "region",
    )
    images = sorted(ts.images or [], key=lambda im: im.section)
    n_raw = len(images)
    if [im.section for im in images] != list(range(n_raw)):
        raise ValueError(f"{run}: tilt image sections are not 0..{n_raw - 1}")
    fr = image_frame(images[0])
    pix = fr.isotropic_spacing
    width, height = fr.size_px
    aligned = {p.z_index for p in hub.per_section_alignment_parameters}
    dark_angles = {im.section: float(im.nominal_tilt_angle or 0.0) for im in images if im.section not in aligned}

    out: Dict[str, Any] = {"run_name": run, "pixel_spacing": pix, "n_raw": n_raw}
    aln_dir = staging / "alignment" / run
    aln_dir.mkdir(parents=True, exist_ok=True)
    if hub.has_x_rotation:
        # IMOD .xf/.tlt/.xtilt (the backend reads that combination; a .aln cannot hold X rotations)
        imod = hub.to_imod(ts_size=(width, height, n_raw), ts_spacing=pix, basename=run)
        for name, obj in (("xf", imod.xf), ("tlt", imod.tlt), ("xtilt", imod.xtilt)):
            p = aln_dir / f"{run}.{name}"
            _refuse_existing(p, overwrite)
            obj.to_file(str(p))
        out["alignment_format"] = "IMOD"
        out["alignment_files"] = [f"{run}.xf", f"{run}.tlt", f"{run}.xtilt"]
        sr.gates.append(
            Gate(
                "x_rotation",
                True,
                value=max(abs(p.volume_x_rotation) for p in hub.per_section_alignment_parameters),
                note="staged as IMOD xf/tlt/xtilt",
            )
        )
        # check: the backend's own reader reproduces the hub
        back = Alignment.from_imod(
            ImodAlignment(
                xf=ImodXF.from_file(aln_dir / f"{run}.xf"),
                tlt=ImodTLT.from_file(aln_dir / f"{run}.tlt"),
                xtilt=ImodXTILT.from_file(aln_dir / f"{run}.xtilt"),
                tiltcom=None,
                newstcom=None,
            )
        )
    else:
        aln = hub.to_aretomo(
            ts_size=(width, height, n_raw),
            dark_angles=dark_angles,
            thickness_px=comp_aln.thickness_px if comp_aln else None,
        )
        aln.header = "# AreTomo Alignment / Priims bprmMn"
        p = aln_dir / f"{run}.aln"
        _refuse_existing(p, overwrite)
        aln.to_file(str(p))
        out["alignment_format"] = "ARETOMO3"
        out["alignment_files"] = [f"{run}.aln"]
        sr.gates.append(Gate("x_rotation", True, value=0.0, note="staged as .aln"))
        back = Alignment.from_aretomo3(AreTomo3ALN.from_file(str(p)))
    worst = _hub_deviation(hub, back)
    sr.gates.append(
        Gate("backend_parser_reproduces_hub", worst < 1e-3, value=worst, expected="< 1e-3 (px / deg; file precision)")
    )

    # rawtlt: nominal angles per raw section
    raw_dir = staging / "rawtlt" / run
    raw_dir.mkdir(parents=True, exist_ok=True)
    rawtlt = raw_dir / f"{run}.rawtlt"
    _refuse_existing(rawtlt, overwrite)
    if any(im.nominal_tilt_angle is None for im in images):
        sr.warnings.append("some tilt images have no nominal_tilt_angle: .rawtlt rows written as 0")
    rawtlt.write_text("".join(f"{float(im.nominal_tilt_angle or 0.0):8.2f}\n" for im in images))
    nominal = [float(im.nominal_tilt_angle or 0.0) for im in images]
    out["tilt_min"], out["tilt_max"] = min(nominal), max(nominal)
    steps = np.diff(sorted(nominal))
    out["tilt_step"] = float(np.round(np.median(steps), 2)) if len(steps) else 0.0
    out["tilt_axis"] = float(hub.get_median_tilt_axis())

    # CTF (CTFFIND-parsable AreTomo3 layout: one header line)
    ctfs = [im.ctf_metadata for im in images]
    if all(c is not None and c.defocus_u is not None for c in ctfs) and not res.optional("no_ctf", absent=False):
        df_hand = res.optional("defocus_hand", companion=(comp_ts.defocus_hand if comp_ts else None))
        rows = [
            cets_ctf.to_aretomo3_row(c, i + 1, df_hand=None if df_hand is None else int(df_hand))
            for i, c in enumerate(ctfs)
        ]
        ctf_dir = staging / "ctf" / run
        ctf_dir.mkdir(parents=True, exist_ok=True)
        p = ctf_dir / f"{run}_CTF.txt"
        _refuse_existing(p, overwrite)
        AreTomo3CTF(rows=rows).to_file(str(p))
        out["ctf"] = True
        sr.gates.append(Gate("ctf_parses_with_one_header_line", _ctffind_parses(p, n_raw), expected=n_raw))
    else:
        out["ctf"] = False
        if any(c is not None for c in ctfs):
            sr.warnings.append("CTF metadata incomplete: no _CTF.txt staged")

    # tilt series / tomogram data links (never promised without data)
    src = _resolve_doc_path(ts.path, doc_dir)
    if src is not None and src.exists() and not str(ts.path).startswith(("http", "s3:")):
        d = staging / "tiltseries" / run
        d.mkdir(parents=True, exist_ok=True)
        _link(src, d / f"{run}{src.suffix}", overwrite)
        out["tiltseries_glob"] = f"tiltseries/{{run_name}}/{{run_name}}{src.suffix}"
    else:
        out["tiltseries_glob"] = None
        out["tiltseries_uri"] = ts.path
    tomo_src = _resolve_doc_path(tomo.path, doc_dir) if tomo.path else None
    out["tomogram"] = None
    if tomo_src is not None and tomo_src.exists() and not str(tomo.path).startswith(("http", "s3:")):
        d = staging / "tomograms" / run
        d.mkdir(parents=True, exist_ok=True)
        _link(tomo_src, d / f"{run}{tomo_src.suffix}", overwrite)
        out["tomogram"] = {
            "glob": f"tomograms/{{run_name}}/{{run_name}}{tomo_src.suffix}",
            "voxel_spacing": reference.spacing_a,
            "size": list(reference.size_px),
            "ctf_corrected": bool(tomo.ctf_corrected),
        }
    else:
        candidates = [
            t
            for t in (region.tomograms or [])
            if t.path
            and not str(t.path).startswith(("http", "s3:"))
            and (_resolve_doc_path(t.path, doc_dir) or Path("/nonexistent")).exists()
        ]
        if candidates:
            t = candidates[0]
            tsrc = _resolve_doc_path(t.path, doc_dir)
            assert tsrc is not None  # candidates were filtered on an existing resolved path
            d = staging / "tomograms" / run
            d.mkdir(parents=True, exist_ok=True)
            _link(tsrc, d / f"{run}{tsrc.suffix}", overwrite)
            tf = image_frame(t)
            out["tomogram"] = {
                "glob": f"tomograms/{{run_name}}/{{run_name}}{tsrc.suffix}",
                "voxel_spacing": tf.isotropic_spacing,
                "size": list(tf.size_px),
                "ctf_corrected": bool(t.ctf_corrected),
                "id": t.id,
            }
        else:
            sr.warnings.append("no local tomogram file to stage: the tomograms block is left as a TODO")
    tcomp = companion.tomograms.get((out["tomogram"] or {}).get("id", tomo.id)) if companion else None
    out["tomogram_meta"] = {
        "processing": tcomp.processing if tcomp else None,
        "reconstruction_method": tcomp.reconstruction_method if tcomp else None,
        "reconstruction_software": tcomp.reconstruction_software if tcomp else None,
    }
    # collection metadata (mdoc) from the companion, when local
    out["mdoc_glob"] = None
    mdoc = (
        _resolve_doc_path(comp_ts.collection_metadata_path, doc_dir)
        if comp_ts and comp_ts.collection_metadata_path
        else None
    )
    if mdoc is not None and mdoc.exists() and not str(comp_ts.collection_metadata_path).startswith(("http", "s3:")):
        d = staging / "collection_metadata" / run
        d.mkdir(parents=True, exist_ok=True)
        _link(mdoc, d / f"{run}.mdoc", overwrite)
        out["mdoc_glob"] = "collection_metadata/{run_name}/*.mdoc"
    elif comp_ts and comp_ts.collection_metadata_path:
        sr.warnings.append(
            f"collection metadata {comp_ts.collection_metadata_path} is not a local file: collection_metadata block left out (the portal requires one per tilt series)"
        )
    else:
        sr.warnings.append(
            "no acquisition mdoc known: collection_metadata block left out (the portal requires one per tilt series)"
        )

    # frames (local movie stacks only)
    frames = []
    if region.movie_stack_collection:
        for series in region.movie_stack_collection.movie_stacks or []:
            for st in series.stacks or []:
                p = _resolve_doc_path(st.path, doc_dir) if st.path else None
                if p is not None and p.exists() and not str(st.path).startswith(("http", "s3:")):
                    frames.append(p)
    if frames:
        d = staging / "frames" / run
        d.mkdir(parents=True, exist_ok=True)
        for p in frames:
            _link(p, d / p.name, overwrite)
        out["frames_glob"] = f"frames/{{run_name}}/*{frames[0].suffix}"
    else:
        out["frames_glob"] = None
    tool = comp_ts.source_tool if comp_ts else None
    out["tilt_alignment_software"] = tool if tool in ("AreTomo3", "Warp") else None
    out["tiltseries_glob"] = res.optional("tiltseries_glob", discovered=out["tiltseries_glob"])
    out["tomograms_glob"] = res.optional("tomograms_glob")
    out["voltage"] = res.optional("voltage", companion=(comp_ts.voltage_kv if comp_ts else None))
    out["cs"] = res.optional("cs", companion=(comp_ts.cs_mm if comp_ts else None))
    out["method_type"] = res.optional("method_type", companion=(comp_aln.method_type if comp_aln else None))
    out["is_portal_standard"] = res.optional(
        "portal_standard", companion=(comp_aln.is_portal_standard if comp_aln else None)
    )
    out["exposure"] = [comp_images[im.id].exposure_dose if im.id in comp_images else None for im in images]
    out["dark_sections"] = sorted(dark_angles)
    sr.outputs.update({"alignment": str(aln_dir), "rawtlt": str(rawtlt)})
    sr.provenance = res.provenance()
    return out


def _refuse_existing(p: Path, overwrite: bool) -> None:
    if p.exists() and not overwrite:
        raise FileExistsError(f"{p} exists (use --overwrite)")


def _link(src: Path, dst: Path, overwrite: bool) -> None:
    if dst.is_symlink() or (dst.exists() and overwrite):
        dst.unlink()
    if not dst.exists():
        os.symlink(src.resolve(), dst)


def _resolve_doc_path(p: Optional[str], doc_dir: Path) -> Optional[Path]:
    if not p:
        return None
    q = Path(p)
    return q if q.is_absolute() else (doc_dir / q)


def _hub_deviation(a: Alignment, b: Alignment) -> float:
    sa, sb = a.sections, b.sections
    if sorted(sa) != sorted(sb):
        return float("inf")
    worst = 0.0
    for z, p in sa.items():
        q = sb[z]
        worst = max(
            worst,
            abs(p.tilt_angle - q.tilt_angle),
            abs(p.tilt_axis_rotation - q.tilt_axis_rotation),
            abs(p.volume_x_rotation - q.volume_x_rotation),
            abs(p.x_offset - q.x_offset),
            abs(p.y_offset - q.y_offset),
        )
    return worst


def _ctffind_parses(path: Path, n_rows: int) -> bool:
    """Emulate the backend's CTFFIND parser: pop exactly one header line, split the rest."""
    lines = path.read_text().strip().splitlines()
    lines.pop(0)
    try:
        rows = [[float(v) for v in ln.split()] for ln in lines]
    except ValueError:
        return False
    return len(rows) == n_rows and all(len(r) in (7, 8) for r in rows)


# ------------------------------------------------------------------------------------------ config draft


def _uniform(runs: List[dict], key: str):
    vals = [r[key] for r in runs]
    return vals[0] if all(v == vals[0] for v in vals) else None


def build_config(
    runs: List[dict], *, deposition_id: int, template: Optional[dict], staging: Path
) -> Tuple[dict, Optional[List[dict]]]:
    """The ingestion-config draft and the run_to_data_map rows (None when every value is uniform)."""
    tpl = template or {}
    tsv_cols: Dict[str, List[Any]] = {}

    def per_run(key: str, kind: str):
        u = _uniform(runs, key)
        if u is not None:
            return u
        tsv_cols[key] = [r[key] for r in runs]
        return f"{kind} {{{key}}}"

    cfg: Dict[str, Any] = {
        "version": CONFIG_VERSION,
        "standardization_config": {
            "deposition_id": int(deposition_id),
            "source_prefix": tpl.get("standardization_config", {}).get("source_prefix", TODO),
        },
        "datasets": tpl.get("datasets")
        or [
            {
                "metadata": {
                    "dataset_identifier": TODO,
                    "dataset_title": TODO,
                    "dataset_description": TODO,
                    "authors": TODO,
                    "dates": TODO,
                    "sample_type": TODO,
                    "organism": TODO,
                    "cross_references": TODO,
                    "funding": TODO,
                    "grid_preparation": TODO,
                    "sample_preparation": TODO,
                },
                "sources": [{"literal": {"value": [TODO]}}],
            }
        ],
        "depositions": [{"sources": [{"literal": {"value": [int(deposition_id)]}}]}],
        "runs": [
            {"sources": [{"source_glob": {"list_glob": "alignment/*", "match_regex": ".*", "name_regex": "(.*)"}}]}
        ],
    }
    ts_meta = copy.deepcopy((tpl.get("tiltseries") or [{}])[0].get("metadata", {}))
    derived = {
        "pixel_spacing": per_run("pixel_spacing", "float"),
        "tilt_axis": per_run("tilt_axis", "float"),
        "tilt_range": {"min": per_run("tilt_min", "float"), "max": per_run("tilt_max", "float")},
        "tilt_step": per_run("tilt_step", "float"),
        "is_aligned": False,
    }
    for k in ("acceleration_voltage", "spherical_aberration_constant"):
        run_key = {"acceleration_voltage": "voltage", "spherical_aberration_constant": "cs"}[k]
        u = _uniform(runs, run_key)
        if u is not None and k not in ts_meta:
            ts_meta[k] = int(round(u * 1000)) if k == "acceleration_voltage" else u
    tool = _uniform(runs, "tilt_alignment_software")
    if tool and "tilt_alignment_software" not in ts_meta:
        ts_meta["tilt_alignment_software"] = tool
    for k in (
        "binning_from_frames",
        "camera",
        "microscope",
        "microscope_optical_setup",
        "data_acquisition_software",
        "tilt_alignment_software",
        "tilt_series_quality",
        "tilting_scheme",
        "total_flux",
    ):
        ts_meta.setdefault(k, TODO)
    ts_meta.update(derived)
    tiltseries_glob = _uniform(runs, "tiltseries_glob")
    cfg["tiltseries"] = [{"metadata": ts_meta, "sources": [{"source_glob": {"list_glob": tiltseries_glob or TODO}}]}]
    frames_meta = copy.deepcopy((tpl.get("frames") or [{}])[0].get("metadata", {}))
    frames_meta.setdefault("dose_rate", TODO)
    frames_meta.setdefault("is_gain_corrected", TODO)
    if _uniform(runs, "frames_glob"):
        cfg["frames"] = [
            {"metadata": frames_meta, "sources": [{"source_glob": {"list_glob": _uniform(runs, "frames_glob")}}]}
        ]
    elif _uniform(runs, "mdoc_glob"):
        # the backend requires a frames block next to collection_metadata; no frame files deposited -> literal default
        cfg["frames"] = [{"metadata": frames_meta, "sources": [{"literal": {"value": ["default"]}}]}]
    if _uniform(runs, "mdoc_glob"):
        cfg["collection_metadata"] = [{"sources": [{"source_glob": {"list_glob": _uniform(runs, "mdoc_glob")}}]}]
    cfg["rawtilts"] = [{"sources": [{"source_glob": {"list_glob": "rawtlt/{run_name}/*.rawtlt"}}]}]
    if all(r["ctf"] for r in runs):
        cfg["ctfs"] = [
            {"metadata": {"format": "CTFFIND"}, "sources": [{"source_glob": {"list_glob": "ctf/{run_name}/*_CTF.txt"}}]}
        ]
    alignments = []
    for fmt in ("ARETOMO3", "IMOD"):
        fmt_runs = [r for r in runs if r["alignment_format"] == fmt]
        if not fmt_runs:
            continue
        method = _uniform(fmt_runs, "method_type") or TODO
        if method not in METHOD_TYPES:
            method = TODO
        standard = _uniform(fmt_runs, "is_portal_standard")
        meta = {
            "format": fmt,
            "alignment_type": "GLOBAL",
            "method_type": method,
            "is_portal_standard": bool(standard) if standard is not None else TODO,
        }
        aln_src: Dict[str, Any]
        if fmt == "ARETOMO3":
            aln_src = {"source_glob": {"list_glob": "alignment/{run_name}/*.aln"}}
        else:
            aln_src = {
                "source_multi_glob": {
                    "list_globs": [
                        "alignment/{run_name}/{run_name}.xf",
                        "alignment/{run_name}/{run_name}.tlt",
                        "alignment/{run_name}/{run_name}.xtilt",
                    ]
                }
            }
        block: Dict[str, Any] = {"metadata": meta, "sources": [aln_src]}
        if len(fmt_runs) != len(runs):
            aln_src["parent_filters"] = {"include": {"run": [f"^{re.escape(r['run_name'])}$" for r in fmt_runs]}}
        alignments.append(block)
    cfg["alignments"] = alignments
    tomo_runs = [r for r in runs if r["tomogram"]]
    if tomo_runs:
        vs = _uniform(tomo_runs, "tomogram")
        voxel = vs["voxel_spacing"] if vs else None
        if voxel is None:
            tsv_cols["voxel_spacing"] = [r["tomogram"]["voxel_spacing"] if r["tomogram"] else "" for r in runs]
        cfg["voxel_spacings"] = [
            {"sources": [{"literal": {"value": [round(voxel, 3)] if voxel else ["float {voxel_spacing}"]}}]}
        ]
        tmeta = copy.deepcopy((tpl.get("tomograms") or [{}])[0].get("metadata", {}))
        meta_src = tomo_runs[0]["tomogram_meta"]
        for k in ("processing", "reconstruction_method", "reconstruction_software"):
            u = _uniform(tomo_runs, "tomogram_meta")
            v = (u or {}).get(k) if u else None
            tmeta.setdefault(k, v if v else TODO)
        tmeta.setdefault("processing_software", tmeta.get("reconstruction_software", TODO))
        tmeta["ctf_corrected"] = _uniform(tomo_runs, "tomogram") and bool(vs["ctf_corrected"]) if vs else TODO
        tmeta.setdefault("align_software", TODO)
        tmeta.setdefault("fiducial_alignment_status", TODO)
        tmeta.setdefault("tomogram_version", 1.0)
        tmeta["voxel_spacing"] = round(voxel, 3) if voxel else "float {voxel_spacing}"
        tmeta.setdefault("offset", {"x": 0, "y": 0, "z": 0})
        tmeta.setdefault("affine_transformation_matrix", [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        tmeta.setdefault("is_visualization_default", True)
        tmeta.setdefault("authors", TODO)
        tmeta.setdefault("dates", TODO)
        del meta_src
        cfg["tomograms"] = [
            {
                "metadata": tmeta,
                "sources": [
                    {
                        "source_glob": {
                            "list_glob": _uniform(runs, "tomograms_glob")
                            or (_uniform(tomo_runs, "tomogram")["glob"] if vs else TODO)
                        }
                    }
                ],
            }
        ]
    else:
        cfg["voxel_spacings"] = [{"sources": [{"literal": {"value": [TODO]}}]}]
        cfg["tomograms"] = TODO
    rows = None
    if tsv_cols:
        rows = [{"run_name": r["run_name"], **{k: tsv_cols[k][i] for k in tsv_cols}} for i, r in enumerate(runs)]
        cfg["standardization_config"]["run_data_map_file"] = "run_to_data_map.tsv"
    return cfg, rows


def write_config(cfg: dict, rows: Optional[List[dict]], staging: Path) -> Path:
    path = staging / "ingestion_config.yaml"
    text = yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True, width=120)
    path.write_text(
        f"# CryoET Data Portal ingestion config DRAFT written by cets-cdp; every '{TODO}' needs curator input\n" + text
    )
    if rows:
        cols = list(rows[0].keys())
        (staging / "run_to_data_map.tsv").write_text(
            "\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)
        )
    return path


def todos(cfg: Any, prefix: str = "") -> List[str]:
    out = []
    if isinstance(cfg, dict):
        for k, v in cfg.items():
            out.extend(todos(v, f"{prefix}.{k}" if prefix else k))
    elif isinstance(cfg, list):
        for i, v in enumerate(cfg):
            out.extend(todos(v, f"{prefix}[{i}]"))
    elif cfg == TODO:
        out.append(prefix)
    return out


# ------------------------------------------------------------------------------------------ validation


def check_sources_resolve(cfg: dict, staging: Path, runs: List[dict], rows: Optional[List[dict]]) -> List[str]:
    """Every glob (with {run_name} and TSV placeholders substituted) must match at least one file per run."""
    problems = []
    table = {r["run_name"]: r for r in (rows or [])}
    for block_name in ("tiltseries", "frames", "collection_metadata", "rawtilts", "ctfs", "alignments", "tomograms"):
        block = cfg.get(block_name)
        if not isinstance(block, list):
            continue
        for entry in block:
            for src in entry.get("sources", []):
                globs = []
                if "source_glob" in src:
                    globs = [src["source_glob"]["list_glob"]]
                elif "source_multi_glob" in src:
                    globs = list(src["source_multi_glob"]["list_globs"])
                include = src.get("parent_filters", {}).get("include", {}).get("run")
                for r in runs:
                    if include and not any(re.match(p, r["run_name"]) for p in include):
                        continue
                    for g in globs:
                        if g == TODO:
                            problems.append(f"{block_name}: source glob is a TODO")
                            continue
                        pat = g.replace("{run_name}", r["run_name"])
                        for k, v in table.get(r["run_name"], {}).items():
                            pat = pat.replace(f"{{{k}}}", str(v))
                        if not list(staging.glob(pat)):
                            problems.append(f"{block_name}: {pat!r} matches nothing for run {r['run_name']}")
    return problems


def _stub_formatted(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _stub_formatted(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_stub_formatted(v) for v in obj]
    if isinstance(obj, str):
        if re.match(r"^float\s*{\w+}\s*$", obj):
            return 1.0
        if re.match(r"^int\s*{\w+}\s*$", obj):
            return 1
    return obj


def check_schema(cfg: dict, backend_path: Optional[str]) -> Tuple[str, List[str]]:
    """In-process validation against the backend's generated pydantic ``Container`` (formatted strings stubbed
    like the backend does — so this proves shape, not that placeholders resolve)."""
    if not backend_path:
        return "skipped", ["CRYOET_DATA_PORTAL_BACKEND_PATH not set"]
    codegen = Path(backend_path) / "schema" / "ingestion_config" / "v1.0.0" / "codegen"
    if not codegen.exists():
        return "skipped", [f"{codegen} not found"]
    sys.path.insert(0, str(codegen))
    try:
        import ingestion_config_models as m  # type: ignore

        m.Container.model_validate(_stub_formatted(cfg))
        return "passed", []
    except Exception as e:  # noqa: BLE001
        return "failed", [str(e)[:4000]]
    finally:
        sys.path.remove(str(codegen))


def check_extended(
    config_path: Path, backend_path: Optional[str], conda_env: Optional[str], out_dir: Path
) -> Tuple[str, List[str]]:
    """The backend's own ``ingestion_config_validate.py`` in its conda env (extended checks)."""
    if not backend_path or not conda_env:
        return "skipped", ["CRYOET_DATA_PORTAL_BACKEND_PATH / CRYOET_DATA_PORTAL_BACKEND_CONDA_ENV not set"]
    script = Path(backend_path) / "schema" / "ingestion_config" / "v1.0.0" / "ingestion_config_validate.py"
    if not script.exists() or shutil.which("conda") is None:
        return "skipped", [f"{script} or conda not available"]
    cmd = ["conda", "run", "-n", conda_env, "python", str(script), str(config_path), "--output-dir", str(out_dir)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600, check=False)
    errors = sorted(out_dir.glob("**/*.json")) if out_dir.exists() else []
    msgs = [(proc.stdout + proc.stderr)[-4000:]] + [str(e) for e in errors]
    return ("passed" if proc.returncode == 0 and not errors else "failed"), msgs


__all__ = [
    "build_config",
    "check_extended",
    "check_schema",
    "check_sources_resolve",
    "stage_run",
    "todos",
    "write_config",
]
