"""Portal run -> CETS ``Region`` (+ companion entries)."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
from cryoet_alignment.io.cets import ctf as cets_ctf
from cryoet_alignment.io.cets.alignment import ReferenceVolume, alignment_to_cets
from cryoet_alignment.io.cets.cli_support import Gate, SeriesReport
from cryoet_alignment.io.cets.companion import (
    AlignmentCompanion,
    AnnotationCompanion,
    ImageCompanion,
    TiltSeriesCompanion,
    TomogramCompanion,
)
from cryoet_alignment.io.cets.config import Resolver
from cryoet_alignment.io.cets.entities import (
    movie_stack_series_entity,
    region_entity,
    tilt_series_entity,
    tomogram_entity,
)
from cryoet_alignment.io.cets.frames import FRAME_CONVENTIONS, image_frame

from cets_cdp.annotations import DEFAULT_SHAPES, portal_annotations_to_cets
from cets_cdp.api import PortalRunData, https_to_s3, implied_voxel, pick_alignment, pick_tomogram


@dataclass
class SeriesResult:
    region: Any
    tilt_series_companion: TiltSeriesCompanion
    alignment_companion: Optional[AlignmentCompanion]
    tomogram_companions: Dict[str, TomogramCompanion]
    annotation_companions: Dict[str, AnnotationCompanion] = field(default_factory=dict)


def _is_identity(m: Optional[list]) -> bool:
    if m is None:
        return True
    a = np.array(m, dtype=np.float64)
    return a.shape == (4, 4) and np.allclose(a, np.eye(4))


def _uri(url: Optional[str], scheme: str) -> Optional[str]:
    if url is None:
        return None
    return https_to_s3(url) if scheme == "s3" else url


def portal_to_cets(
    data: PortalRunData,
    res: Resolver,
    sr: SeriesReport,
    *,
    alignment_id: Optional[int] = None,
    voxel: Optional[float] = None,
    annotations: bool = True,
    annotation_shapes: Sequence[str] = DEFAULT_SHAPES,
) -> SeriesResult:
    stem = data.run_name
    scheme = res.value("uri_scheme", default="https")
    pix = data.pixel_spacing
    res.resolve("pix", discovered=pix, note="TiltSeries.pixel_spacing")
    width, height, n_raw = data.size
    if data.is_aligned:
        raise ValueError(
            f"{stem}: the portal tilt series is flagged is_aligned - the per-section shifts do not apply to an aligned stack"
        )

    sections = {s.z_index: s for s in data.sections}
    if sorted(sections) != list(range(n_raw)):
        sr.warnings.append(f"per-section parameters cover {sorted(sections)} of {n_raw} raw sections")
    nominal = [sections[z].raw_angle if z in sections else 0.0 for z in range(n_raw)]
    if any(z not in sections for z in range(n_raw)):
        sr.warnings.append("sections without PerSectionParameters get nominal_tilt_angle 0 (no source)")
    doses = [sections[z].accumulated_dose if z in sections else None for z in range(n_raw)]
    ctfs = None
    if all(z in sections and sections[z].major_defocus_a is not None for z in range(n_raw)):
        ctfs = [
            cets_ctf.from_portal_values(
                sections[z].major_defocus_a,
                sections[z].minor_defocus_a,
                sections[z].astigmatic_angle_deg,
                sections[z].phase_shift_rad,
            )
            for z in range(n_raw)
        ]
    frame_urls = [sections[z].frame_url if z in sections else None for z in range(n_raw)]
    have_frames = all(frame_urls)
    movie_ids = [f"{stem}_movie_{z}" for z in range(n_raw)] if have_frames else None

    ts = tilt_series_entity(
        tilt_series_id=stem,
        path=_uri(data.https_mrc_file, scheme),
        width=width,
        height=height,
        pixel_size_a=pix,
        nominal_angles=nominal,
        doses=doses,
        ctfs=ctfs,
        movie_stack_ids=movie_ids,
        movie_stack_series_id=f"{stem}_movies" if have_frames else None,
    )
    movie_series = []
    if have_frames and movie_ids is not None:
        movie_series.append(
            movie_stack_series_entity(
                series_id=f"{stem}_movies",
                stacks=[{"id": movie_ids[z], "path": _uri(frame_urls[z], scheme)} for z in range(n_raw)],
            )
        )

    # tomograms: every portal tomogram at the voxel spacing the portal declares (the portal's own coordinate
    # convention; the value implied by the raw field of view is recorded in the companion for information only)
    ref_portal = pick_tomogram(data, voxel)
    tomograms: List[Any] = []
    tomo_comps: Dict[str, TomogramCompanion] = {}
    ref_tomo = None
    tomograms_by_portal_id: Dict[int, Any] = {}
    for t in data.tomograms:
        tomo = tomogram_entity(
            tomogram_id=f"{stem}_tomo_{t.id}",
            path=_uri(t.https_mrc_file or t.https_omezarr_dir, scheme),
            size_px=t.size,
            voxel_size_a=t.voxel_spacing,
            tilt_series_id=stem,
            ctf_corrected=bool(t.ctf_corrected),
        )
        tomograms.append(tomo)
        tomograms_by_portal_id[t.id] = tomo
        tomo_comps[tomo.id] = TomogramCompanion(
            voxel_header_a=t.voxel_spacing,
            voxel_implied_a=implied_voxel(data, t),
            processing=t.processing,
            reconstruction_method=t.reconstruction_method,
            reconstruction_software=t.reconstruction_software,
            source_ref=f"portal tomogram {t.id}",
        )
        if ref_portal is not None and t.id == ref_portal.id:
            ref_tomo = tomo
    aln = pick_alignment(data, alignment_id)
    if ref_tomo is not None and ref_portal is not None:
        res.resolve(
            "reference_tomogram",
            discovered=ref_tomo.id,
            note=f"portal tomogram {ref_portal.id} at {ref_portal.voxel_spacing} Å"
            + (" (--voxel)" if voxel else " (portal standard / finest)"),
        )
    elif aln is not None and all(aln.volume_dimension_a.get(k) for k in "xyz"):
        # no reconstruction on the portal: the alignment's own box (portal volume_dimension) at the tilt pixel
        box = aln.volume_dimension_a
        size = tuple(int(round(float(box[k]) / pix)) for k in "xyz")
        ref_tomo = tomogram_entity(
            tomogram_id=f"{stem}_volume", path=None, size_px=size, voxel_size_a=pix, tilt_series_id=stem
        )
        tomograms.insert(0, ref_tomo)
        tomo_comps[ref_tomo.id] = TomogramCompanion(
            voxel_implied_a=pix, source_ref=f"portal alignment {aln.id} volume_dimension at pixel_spacing; no file"
        )
        res.resolve(
            "reference_tomogram", discovered=ref_tomo.id, note="no tomogram on the portal: alignment volume_dimension"
        )
        sr.warnings.append(
            "no tomogram on the portal for this run: the alignment's volume_dimension box is the reference frame"
        )

    alignments = []
    aln_comp = None
    if aln is not None and ref_tomo is not None:
        hub = aln.hub
        if hub is None:
            raise ValueError(
                f"{stem}: portal alignment {aln.id} has no alignment_metadata.json (cached record without it)"
            )
        if not _is_identity(aln.affine_transformation_matrix) or any(
            abs(float(v or 0.0)) > 0 for v in aln.volume_offset_a.values()
        ):
            raise ValueError(
                f"{stem}: portal alignment {aln.id} carries a non-identity affine_transformation_matrix or a volume_offset; "
                "these registration terms are not modelled by the rigid profile",
            )
        if aln.tilt_offset:
            sr.warnings.append(
                f"alignment tilt_offset {aln.tilt_offset} deg is already inside the per-section tilt angles (AreTomo3 semantics); recorded only"
            )
        ref_frame = image_frame(ref_tomo)
        native = {k: float(n * s) for k, n, s in zip("xyz", ref_frame.size_px, ref_frame.spacing_a, strict=True)}
        portal_box = {k: float(v) for k, v in aln.volume_dimension_a.items() if v is not None}
        if portal_box:
            dev = max(abs(portal_box.get(k, native[k]) - native[k]) for k in "xyz")
            sr.gates.append(
                Gate(
                    "portal_volume_box_vs_reference_tomogram",
                    dev <= ref_frame.isotropic_spacing,
                    value=portal_box,
                    expected=native,
                    note=f"deviation {dev:.2f} Å; the portal computes its box from the ingested tomogram's rounded voxel",
                )
            )
        hub.volume_dimension = native
        cets_alignment = alignment_to_cets(
            hub,
            tilt_series_id=stem,
            alignment_name=f"portal{aln.id}",
            image=image_frame(ts.images[0]),
            reference=ReferenceVolume.from_tomogram(ref_tomo),
            frame=FRAME_CONVENTIONS["ARETOMO3"],
        )
        alignments.append(cets_alignment)
        dropped = []
        if (aln.alignment_type or "").upper() == "LOCAL":
            dropped.append("LOCAL alignment: the portal metadata carries the rigid per-section parameters only")
        sr.dropped.extend(dropped)
        sr.gates.append(
            Gate(
                "rows",
                len(cets_alignment.projection_alignments) == len(hub.per_section_alignment_parameters),
                value=len(cets_alignment.projection_alignments),
            )
        )
        aln_comp = AlignmentCompanion(
            name=f"portal{aln.id}",
            tilt_series_id=stem,
            format=hub.format,
            alignment_type=aln.alignment_type,
            method_type=aln.alignment_method,
            is_portal_standard=aln.is_portal_standard,
            reference_tomogram_id=ref_tomo.id,
            tomogram_ids=[t.id for t in tomograms],
            native_volume_dimension_a=native,
            frame_convention={"image_center": "half", "volume_center": "half"},
            dropped=dropped,
            source_ref=f"portal:{data.dataset_id}/{stem} alignment {aln.id}; portal volume_dimension {portal_box}",
        )
    elif aln is None:
        sr.warnings.append("no alignment on the portal for this run")
    else:
        sr.warnings.append("no tomogram and no alignment box on the portal: alignment not encoded")

    images = {}
    for z in range(n_raw):
        s = sections.get(z)
        images[ts.images[z].id] = ImageCompanion(
            acquisition_index_1b=s.acquisition_order_1b if s else None,
            exposure_dose=s.exposure_dose if s else None,
            stage_angle_deg=s.raw_angle if s else None,
            frame_name=(s.frame_url.rsplit("/", 1)[-1] if s and s.frame_url else None),
            use_tilt=(z in aln.hub.sections) if aln is not None and aln.hub is not None else None,
        )
    ts_comp = TiltSeriesCompanion(
        source_tool="cryoET Data Portal",
        voltage_kv=res.optional("voltage", discovered=data.voltage_kv, note="TiltSeries.acceleration_voltage"),
        cs_mm=res.optional("cs", discovered=data.cs_mm, note="TiltSeries.spherical_aberration_constant"),
        amplitude_contrast=res.optional("amp_contrast"),
        pixel_size_acquisition_a=pix,
        tilt_axis_nominal_deg=data.tilt_axis_deg,
        alpha_offset_deg=aln.tilt_offset if aln else None,
        beta_offset_deg=aln.x_rotation_offset if aln else None,
        collection_metadata_path=_uri(data.mdoc_url, scheme),
        images=images,
    )
    ann_entities: List[Any] = []
    ann_comps: Dict[str, AnnotationCompanion] = {}
    if annotations and data.annotations:
        imported = portal_annotations_to_cets(
            data,
            tomograms_by_portal_id,
            sr,
            alignment_id=aln.id if aln is not None else None,
            shapes=annotation_shapes,
            scheme=scheme,
        )
        ann_entities = imported.entities
        ann_comps = imported.companions
        res.resolve(
            "annotations",
            discovered=len(ann_entities),
            note=f"of {sum(len(a.files) for a in data.annotations)} portal annotation files",
        )
    region = region_entity(
        region_id=stem,
        tilt_series=[ts],
        alignments=alignments,
        tomograms=tomograms,
        movie_stack_series=movie_series,
        annotations=ann_entities,
    )
    sr.provenance = res.provenance()
    return SeriesResult(region, ts_comp, aln_comp, tomo_comps, ann_comps)
