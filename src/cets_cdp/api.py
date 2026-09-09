"""Portal API access: everything one run needs, as plain data (buildable without a client for tests).

Fields (client 4.8.0, verified on 10445 / TS_105_5): ``TiltSeries{pixel_spacing, size_x/y/z, tilt_axis,
acceleration_voltage (V), spherical_aberration_constant (mm), is_aligned, https_mrc_file, https_omezarr_dir}``,
``Frame{acquisition_order (0-based), exposure_dose, accumulated_dose (before the frame), https_frame_path}``,
``PerSectionParameters{z_index, raw_angle, frame_id, major/minor_defocus (Å), astigmatic_angle (deg),
phase_shift (RADIANS)}``, ``Alignment{id, alignment_type, alignment_method, is_portal_standard,
volume_*_dimension/offset (Å), x_rotation_offset, tilt_offset, affine_transformation_matrix (string),
https_alignment_metadata}`` and ``Tomogram{id, voxel_spacing, size_x/y/z, processing, ctf_corrected,
reconstruction_method, reconstruction_software, is_portal_standard, https_mrc_file, https_omezarr_dir}``.
The hub alignment is built from the ``alignment_metadata.json`` the ingestion wrote (identical to what the
API serves; per-section offsets are in PIXELS despite the schema docstring).
"""

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import requests
from cryoet_alignment.io.cryoet_data_portal import Alignment

FILES_HOST = "https://files.cryoetdataportal.cziscience.com/"
S3_BUCKET = "s3://cryoet-data-portal-public/"

PORTAL_SOURCE_FORMS = "portal:<dataset-id>[/<run-name>][@alignment=ID,voxel=A]"
_PORTAL_RE = re.compile(r"^portal:(?P<dataset>\d+)(?:/(?P<run>[^@]+))?(?:@(?P<opts>.*))?$")


@dataclass(frozen=True)
class PortalSource:
    dataset_id: int
    run_name: Optional[str] = None
    alignment_id: Optional[int] = None
    voxel_spacing: Optional[float] = None


def parse_portal_source(token: str) -> PortalSource:
    m = _PORTAL_RE.match(token.strip())
    if not m:
        raise ValueError(f"{token!r}: portal sources are {PORTAL_SOURCE_FORMS}")
    alignment = voxel = None
    for part in (m.group("opts") or "").split(","):
        part = part.strip()
        if not part:
            continue
        key, _, val = part.partition("=")
        if key == "alignment" and val:
            alignment = int(val)
        elif key == "voxel" and val:
            voxel = float(val)
        else:
            raise ValueError(f"{token!r}: unknown portal option {part!r} (alignment=ID, voxel=A)")
    run = m.group("run")
    return PortalSource(int(m.group("dataset")), run.strip() if run else None, alignment, voxel)


def https_to_s3(url: Optional[str]) -> Optional[str]:
    return S3_BUCKET + url[len(FILES_HOST) :] if url and url.startswith(FILES_HOST) else url


@dataclass
class PortalSection:
    z_index: int
    raw_angle: float
    acquisition_order_1b: Optional[int]
    exposure_dose: Optional[float]
    accumulated_dose: Optional[float]
    frame_url: Optional[str]
    major_defocus_a: Optional[float] = None
    minor_defocus_a: Optional[float] = None
    astigmatic_angle_deg: Optional[float] = None
    phase_shift_rad: Optional[float] = None


@dataclass
class PortalTomogram:
    id: int
    voxel_spacing: float
    size: tuple
    processing: Optional[str]
    ctf_corrected: Optional[bool]
    reconstruction_method: Optional[str]
    reconstruction_software: Optional[str]
    is_portal_standard: Optional[bool]
    https_mrc_file: Optional[str]
    https_omezarr_dir: Optional[str]


@dataclass
class PortalAlignment:
    id: int
    alignment_type: Optional[str]
    alignment_method: Optional[str]
    is_portal_standard: Optional[bool]
    https_alignment_metadata: str
    volume_dimension_a: dict
    volume_offset_a: dict
    tilt_offset: float
    x_rotation_offset: float
    affine_transformation_matrix: Optional[list]
    hub: Optional[Alignment] = None
    metadata: Optional[dict] = None


@dataclass
class PortalRunData:
    dataset_id: int
    run_id: int
    run_name: str
    tiltseries_id: int
    pixel_spacing: float
    size: tuple
    is_aligned: bool
    voltage_kv: Optional[float]
    cs_mm: Optional[float]
    tilt_axis_deg: Optional[float]
    https_mrc_file: Optional[str]
    https_omezarr_dir: Optional[str]
    mdoc_url: Optional[str]
    sections: List[PortalSection] = field(default_factory=list)
    alignments: List[PortalAlignment] = field(default_factory=list)
    tomograms: List[PortalTomogram] = field(default_factory=list)

    @property
    def stem(self) -> str:
        return self.run_name


def portal_client(url: Optional[str] = None):
    from cryoet_data_portal import Client

    return Client(url) if url else Client()


def find_runs(client, source: PortalSource) -> list:
    from cryoet_data_portal import Run

    if source.run_name:
        hits = Run.find(client, [Run.dataset_id == source.dataset_id, Run.name == source.run_name])
        if not hits:
            raise ValueError(f"portal run {source.run_name!r} not found in dataset {source.dataset_id}")
        return hits
    return Run.find(client, [Run.dataset_id == source.dataset_id])


def fetch_alignment_metadata(url: str, cache: Optional[Path] = None) -> dict:
    if cache is not None and cache.exists():
        return json.loads(cache.read_text())
    j = requests.get(url, timeout=120).json()
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(j, indent=1))
    return j


def fetch_run(client, run, *, cache_dir: Optional[Path] = None) -> PortalRunData:
    from cryoet_data_portal import Alignment as ApiAlignment
    from cryoet_data_portal import Frame, FrameAcquisitionFile, PerSectionParameters, TiltSeries, Tomogram

    ts_list = TiltSeries.find(client, [TiltSeries.run_id == run.id])
    if not ts_list:
        raise ValueError(f"portal run {run.name} ({run.id}) has no tilt series")
    ts = ts_list[0]
    frames = {f.id: f for f in Frame.find(client, [Frame.run_id == run.id])}
    psp = sorted(PerSectionParameters.find(client, [PerSectionParameters.run_id == run.id]), key=lambda p: p.z_index)
    sections = []
    for p in psp:
        f = frames.get(p.frame_id)
        sections.append(
            PortalSection(
                z_index=int(p.z_index),
                raw_angle=float(p.raw_angle),
                acquisition_order_1b=(int(f.acquisition_order) + 1)
                if f is not None and f.acquisition_order is not None
                else None,
                exposure_dose=float(f.exposure_dose) if f is not None and f.exposure_dose is not None else None,
                accumulated_dose=float(f.accumulated_dose)
                if f is not None and f.accumulated_dose is not None
                else None,
                frame_url=f.https_frame_path if f is not None else None,
                major_defocus_a=p.major_defocus,
                minor_defocus_a=p.minor_defocus,
                astigmatic_angle_deg=p.astigmatic_angle,
                phase_shift_rad=p.phase_shift,
            )
        )
    alignments = []
    for a in ApiAlignment.find(client, [ApiAlignment.run_id == run.id]):
        cache = (cache_dir / run.name / f"alignment_{a.id}.json") if cache_dir else None
        meta = fetch_alignment_metadata(a.https_alignment_metadata, cache)
        affine = None
        if a.affine_transformation_matrix:
            try:
                affine = (
                    json.loads(a.affine_transformation_matrix)
                    if isinstance(a.affine_transformation_matrix, str)
                    else a.affine_transformation_matrix
                )
            except ValueError:
                affine = None
        alignments.append(
            PortalAlignment(
                id=int(a.id),
                alignment_type=a.alignment_type,
                alignment_method=a.alignment_method,
                is_portal_standard=getattr(a, "is_portal_standard", None),
                https_alignment_metadata=a.https_alignment_metadata,
                volume_dimension_a={"x": a.volume_x_dimension, "y": a.volume_y_dimension, "z": a.volume_z_dimension},
                volume_offset_a={"x": a.volume_x_offset, "y": a.volume_y_offset, "z": a.volume_z_offset},
                tilt_offset=float(a.tilt_offset or 0.0),
                x_rotation_offset=float(a.x_rotation_offset or 0.0),
                affine_transformation_matrix=affine or meta.get("affine_transformation_matrix"),
                hub=Alignment(**{k: meta[k] for k in Alignment.model_fields if k in meta}),
                metadata=meta,
            )
        )
    mdocs = FrameAcquisitionFile.find(client, [FrameAcquisitionFile.run_id == run.id])
    tomos = [
        PortalTomogram(
            id=int(t.id),
            voxel_spacing=float(t.voxel_spacing),
            size=(int(t.size_x), int(t.size_y), int(t.size_z)),
            processing=t.processing,
            ctf_corrected=t.ctf_corrected,
            reconstruction_method=t.reconstruction_method,
            reconstruction_software=t.reconstruction_software,
            is_portal_standard=getattr(t, "is_portal_standard", None),
            https_mrc_file=t.https_mrc_file,
            https_omezarr_dir=t.https_omezarr_dir,
        )
        for t in Tomogram.find(client, [Tomogram.run_id == run.id])
    ]
    return PortalRunData(
        dataset_id=int(run.dataset_id),
        run_id=int(run.id),
        run_name=str(run.name),
        tiltseries_id=int(ts.id),
        pixel_spacing=float(ts.pixel_spacing),
        size=(int(ts.size_x), int(ts.size_y), int(ts.size_z)),
        is_aligned=bool(ts.is_aligned),
        voltage_kv=float(ts.acceleration_voltage) / 1000.0 if ts.acceleration_voltage else None,
        cs_mm=float(ts.spherical_aberration_constant) if ts.spherical_aberration_constant is not None else None,
        tilt_axis_deg=float(ts.tilt_axis) if ts.tilt_axis is not None else None,
        https_mrc_file=ts.https_mrc_file,
        https_omezarr_dir=ts.https_omezarr_dir,
        mdoc_url=mdocs[0].https_mdoc_path if mdocs else None,
        sections=sections,
        alignments=alignments,
        tomograms=tomos,
    )


def implied_voxel(data: PortalRunData, tomo: PortalTomogram) -> float:
    """Raw-field-implied voxel of a portal tomogram (its header voxel is rounded to 2-3 decimals)."""
    return data.pixel_spacing * data.size[0] / tomo.size[0]


def pick_tomogram(data: PortalRunData, voxel_spacing: Optional[float] = None) -> Optional[PortalTomogram]:
    if not data.tomograms:
        return None
    if voxel_spacing is not None:
        cands = [t for t in data.tomograms if abs(t.voxel_spacing - voxel_spacing) < 1e-3]
        if not cands:
            raise ValueError(
                f"no portal tomogram at voxel spacing {voxel_spacing} (has {sorted({t.voxel_spacing for t in data.tomograms})})"
            )
        return cands[0]
    # default: the portal-standard tomogram, finest voxel first, lowest id last
    return sorted(data.tomograms, key=lambda t: (not bool(t.is_portal_standard), t.voxel_spacing, t.id))[0]


def pick_alignment(data: PortalRunData, alignment_id: Optional[int] = None) -> Optional[PortalAlignment]:
    if not data.alignments:
        return None
    if alignment_id is not None:
        for a in data.alignments:
            if a.id == alignment_id:
                return a
        raise ValueError(f"alignment {alignment_id} not in run {data.run_name} (has {[a.id for a in data.alignments]})")
    ranked = sorted(data.alignments, key=lambda a: (not bool(a.is_portal_standard), a.id))
    return ranked[0]


# ------------------------------------------------------------------ recording / replay (tests, offline use)


def run_data_to_json(data: PortalRunData) -> dict:
    from dataclasses import asdict

    d = asdict(data)
    for a, src in zip(d["alignments"], data.alignments, strict=True):
        a["hub"] = json.loads(str(src.hub)) if src.hub is not None else None
    return d


def run_data_from_json(d: dict) -> PortalRunData:
    d = dict(d)
    sections = [PortalSection(**s) for s in d.pop("sections")]
    alignments = []
    for a in d.pop("alignments"):
        a = dict(a)
        hub = a.pop("hub", None)
        a["hub"] = Alignment(**hub) if hub else None
        alignments.append(PortalAlignment(**a))
    tomograms = [PortalTomogram(**{**t, "size": tuple(t["size"])}) for t in d.pop("tomograms")]
    d["size"] = tuple(d["size"])
    return PortalRunData(**d, sections=sections, alignments=alignments, tomograms=tomograms)
