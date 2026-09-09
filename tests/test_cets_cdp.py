"""cets-cryoet-data-portal gates on the recorded 10445 / TS_105_5 run (alignment 18924) plus an optional live
run (``CETS_CDP_NETWORK=1``) and the backend validators (``CRYOET_DATA_PORTAL_BACKEND_PATH``)."""

import json
import os
from pathlib import Path

import numpy as np
import pytest
import yaml
from click.testing import CliRunner
from cryoet_alignment.io.aretomo3 import AreTomo3ALN
from cryoet_alignment.io.cets.cli_support import SeriesReport
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.config import Resolver
from cryoet_alignment.io.cets.entities import dataset_entity, dump_json, load_dataset, validate_document

from cets_cdp.annotations import read_ndjson_points, relion4_star_to_points
from cets_cdp.api import PortalRunData, parse_portal_source, run_data_from_json
from cets_cdp.cli import main
from cets_cdp.from_cets import TODO, build_config, check_schema, check_sources_resolve, todos
from cets_cdp.to_cets import portal_to_cets

DATA = Path(__file__).parent / "data"
RUN = "TS_105_5"


def _run(args, expect_ok=True):
    result = CliRunner().invoke(main, args, catch_exceptions=False)
    if expect_ok:
        assert result.exit_code == 0, result.stdout + result.stderr
    return result


@pytest.fixture
def data() -> PortalRunData:
    d = run_data_from_json(json.loads((DATA / "10445_TS_105_5.json").read_text()))
    for a in d.annotations:  # recorded ndjson files live under tests/data
        for f in a.files:
            if f.local_path and not Path(f.local_path).is_absolute():
                f.local_path = str(DATA / f.local_path)
    return d


def _to_cets(data: PortalRunData, out: Path, **kw):
    sr = SeriesReport(RUN)
    res = Resolver("cets-cdp", "to-cets", cli={}, warn=sr.warnings.append)
    result = portal_to_cets(data, res, sr, **kw)
    ds = dataset_entity("10445", [result.region])
    validate_document(ds)
    dump_json(ds, out)
    comp = Companion(
        generator="test",
        tilt_series={RUN: result.tilt_series_companion},
        alignments=[result.alignment_companion],
        tomograms=result.tomogram_companions,
        annotations=result.annotation_companions,
    )
    comp.dump(Companion.path_for(out))
    return ds, comp, sr


def test_source_spec():
    s = parse_portal_source("portal:10445/TS_105_5@alignment=18924,voxel=4.99")
    assert (s.dataset_id, s.run_name, s.alignment_id, s.voxel_spacing) == (10445, "TS_105_5", 18924, 4.99)
    assert parse_portal_source("portal:10445").run_name is None
    with pytest.raises(ValueError):
        parse_portal_source("portal:x")


def test_recorded_run_to_cets(data, tmp_path):
    ds, comp, sr = _to_cets(data, tmp_path / "10445.cets.json")
    region = ds.regions[0]
    ts = region.tilt_series[0]
    assert len(ts.images) == 31 and ts.images[0].width == 4096
    assert ts.images[0].nominal_tilt_angle == pytest.approx(data.sections[0].raw_angle) and ts.images[
        0
    ].ctf_metadata.defocus_u == pytest.approx(data.sections[0].major_defocus_a)
    assert ts.images[0].ctf_metadata.phase_shift == pytest.approx(np.degrees(data.sections[0].phase_shift_rad or 0.0))
    # tomograms carry the portal's declared voxel spacing; the reference is the portal-standard one (finest)
    ids = [t.id for t in region.tomograms]
    assert len(ids) == 4 and not any(i.endswith("_volume") for i in ids)
    for t, tc in zip(region.tomograms, [comp.tomograms[i] for i in ids], strict=True):
        assert t.coordinate_transformations[0].sequence[1].scale[0] == tc.voxel_header_a
        assert tc.voxel_header_a in (4.99, 10.012) and tc.processing
        assert tc.voxel_implied_a == pytest.approx(4096 * 1.54 / t.width)  # information only
    ref = comp.alignments[0].reference_tomogram_id
    ref_tomo = next(t for t in region.tomograms if t.id == ref)
    assert ref_tomo.coordinate_transformations[0].sequence[1].scale[0] == 4.99 and ref_tomo.width == 1260
    assert comp.alignments[0].alignment_type == "LOCAL" and "rigid per-section" in comp.alignments[0].dropped[0]
    # the alignment box equals the portal's own volume_dimension (size x declared voxel)
    assert comp.alignments[0].method_type == "projection_matching"
    assert comp.alignments[0].native_volume_dimension_a == pytest.approx(
        {"x": 1260 * 4.99, "y": 1260 * 4.99, "z": 368 * 4.99}
    )
    assert data.alignments[0].volume_dimension_a["x"] == pytest.approx(1260 * 4.99)
    assert comp.tilt_series[RUN].voltage_kv == 300.0 and comp.tilt_series[RUN].tilt_axis_nominal_deg == -96.0
    assert comp.tilt_series[RUN].images[f"{RUN}_0"].acquisition_index_1b is not None
    pa = region.alignments[0].projection_alignments[0]
    assert pa.id == f"{RUN}_portal18924_align_0" and pa.sequence[2].translation == pytest.approx(
        [-16.81 * 1.54, -79.316 * 1.54]
    )
    assert any(g.name == "portal_volume_box_vs_reference_tomogram" and g.passed for g in sr.gates)


def test_roundtrip_to_aln_equals_portal_aln(data, tmp_path):
    out = tmp_path / "10445.cets.json"
    _to_cets(data, out)
    r = _run(["from-cets", str(out), "-o", str(tmp_path / "stage"), "--deposition-id", "1", "--no-validate"])
    assert "[ok ] backend_parser_reproduces_hub" in r.stdout
    got = AreTomo3ALN.from_file(tmp_path / "stage" / "alignment" / RUN / f"{RUN}.aln")
    portal = AreTomo3ALN.from_file(DATA / "10445_TS_105_5.aln")
    assert [str(g) for g in got.GlobalAlignments] == [str(g) for g in portal.GlobalAlignments]
    assert (tmp_path / "stage" / "ctf" / RUN / f"{RUN}_CTF.txt").exists()
    cfg = yaml.safe_load((tmp_path / "stage" / "ingestion_config.yaml").read_text())
    assert cfg["alignments"][0]["metadata"] == {
        "format": "ARETOMO3",
        "alignment_type": "GLOBAL",
        "method_type": "projection_matching",
        "is_portal_standard": True,
    }
    ts_meta = cfg["tiltseries"][0]["metadata"]
    assert (
        ts_meta["tilt_axis"] == pytest.approx(-96.3299)
        and ts_meta["tilt_range"] == {"min": -45.03, "max": 44.96}
        and ts_meta["tilt_step"] == 3.0
    )
    assert ts_meta["pixel_spacing"] == 1.54 and ts_meta["acceleration_voltage"] == 300000
    assert ts_meta["camera"] == TODO and cfg["tomograms"] == TODO  # nothing inferred, nothing promised without data
    assert "collection_metadata" not in cfg  # the mdoc is an https URL, not local
    assert "TODO(curator): tomograms" in r.stderr


def test_refusals(data, tmp_path):
    aligned = run_data_from_json({**json.loads((DATA / "10445_TS_105_5.json").read_text()), "is_aligned": True})
    with pytest.raises(ValueError, match="is_aligned"):
        _to_cets(aligned, tmp_path / "a.cets.json")
    d2 = run_data_from_json(json.loads((DATA / "10445_TS_105_5.json").read_text()))
    d2.alignments[0].volume_offset_a = {"x": 0.0, "y": 0.0, "z": 12.0}
    with pytest.raises(ValueError, match="volume_offset"):
        _to_cets(d2, tmp_path / "b.cets.json")
    d3 = run_data_from_json(json.loads((DATA / "10445_TS_105_5.json").read_text()))
    d3.alignments[0].affine_transformation_matrix = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 625, 1]]
    with pytest.raises(ValueError, match="affine_transformation_matrix"):
        _to_cets(d3, tmp_path / "c.cets.json")


def test_x_rotation_stages_imod(data, tmp_path):
    from cryoet_alignment.io.cets.rotation import tilt_matrix

    out = tmp_path / "10445.cets.json"
    _to_cets(data, out)
    doc = json.loads(out.read_text())
    for pa in doc["regions"][0]["alignments"][0]["projection_alignments"]:
        tilt = float(np.degrees(np.arcsin(-pa["sequence"][0]["affine"][2][0])))
        pa["sequence"][0]["affine"] = tilt_matrix(tilt, 0.4).tolist()
    out.write_text(json.dumps(doc))
    r = _run(["from-cets", str(out), "-o", str(tmp_path / "stage"), "--deposition-id", "1", "--no-validate"])
    assert "staged as IMOD xf/tlt/xtilt" in r.stdout and "[ok ] backend_parser_reproduces_hub" in r.stdout
    assert sorted(p.name for p in (tmp_path / "stage" / "alignment" / RUN).iterdir()) == [
        f"{RUN}.tlt",
        f"{RUN}.xf",
        f"{RUN}.xtilt",
    ]
    cfg = yaml.safe_load((tmp_path / "stage" / "ingestion_config.yaml").read_text())
    assert (
        cfg["alignments"][0]["metadata"]["format"] == "IMOD"
        and "source_multi_glob" in cfg["alignments"][0]["sources"][0]
    )
    xtilt = [float(v) for v in (tmp_path / "stage" / "alignment" / RUN / f"{RUN}.xtilt").read_text().split()]
    assert xtilt == pytest.approx([0.4] * 31, abs=1e-6)


def test_config_heterogeneous_runs_use_run_data_map(tmp_path):
    runs = []
    for name, pix, axis in (("A", 1.54, -96.0), ("B", 1.34, -95.0)):
        runs.append(
            {
                "run_name": name,
                "pixel_spacing": pix,
                "n_raw": 31,
                "alignment_format": "ARETOMO3",
                "alignment_files": [f"{name}.aln"],
                "tilt_min": -45.0,
                "tilt_max": 45.0,
                "tilt_step": 3.0,
                "tilt_axis": axis,
                "ctf": True,
                "tiltseries_glob": "tiltseries/{run_name}/{run_name}.mrc",
                "tomogram": {
                    "glob": "tomograms/{run_name}/{run_name}.mrc",
                    "voxel_spacing": 6.16,
                    "size": [1024, 1024, 500],
                    "ctf_corrected": False,
                },
                "tomogram_meta": {
                    "processing": "raw",
                    "reconstruction_method": "WBP",
                    "reconstruction_software": "AreTomo3",
                },
                "frames_glob": None,
                "mdoc_glob": None,
                "voltage": 300.0,
                "cs": 2.7,
                "method_type": "projection_matching",
                "is_portal_standard": True,
                "exposure": [None] * 31,
                "dark_sections": [],
                "tilt_alignment_software": "AreTomo3",
                "tomograms_glob": None,
            }
        )
    cfg, rows = build_config(runs, deposition_id=7, template=None, staging=tmp_path)
    assert rows == [
        {"run_name": "A", "pixel_spacing": 1.54, "tilt_axis": -96.0},
        {"run_name": "B", "pixel_spacing": 1.34, "tilt_axis": -95.0},
    ]
    assert (
        cfg["tiltseries"][0]["metadata"]["pixel_spacing"] == "float {pixel_spacing}"
        and cfg["standardization_config"]["run_data_map_file"] == "run_to_data_map.tsv"
    )
    assert cfg["tomograms"][0]["metadata"]["processing"] == "raw" and cfg["voxel_spacings"][0]["sources"][0]["literal"][
        "value"
    ] == [6.16]
    assert "tiltseries[0].metadata.camera" in todos(cfg)
    assert check_sources_resolve(cfg, tmp_path, runs, rows)  # nothing staged here -> problems reported


@pytest.mark.skipif(not os.environ.get("CRYOET_DATA_PORTAL_BACKEND_PATH"), reason="backend checkout not configured")
def test_schema_check_with_backend(tmp_path):
    cfg = {
        "version": "1.1.0",
        "standardization_config": {"deposition_id": 1, "source_prefix": "x"},
        "datasets": [{"metadata": {"dataset_identifier": 1}, "sources": [{"literal": {"value": [1]}}]}],
        "depositions": [{"sources": [{"literal": {"value": [1]}}]}],
        "runs": [{"sources": [{"literal": {"value": ["A"]}}]}],
        "voxel_spacings": [{"sources": [{"literal": {"value": [6.16]}}]}],
    }
    status, msgs = check_schema(cfg, os.environ["CRYOET_DATA_PORTAL_BACKEND_PATH"])
    assert status in ("passed", "failed") and isinstance(msgs, list)


@pytest.mark.skipif(not os.environ.get("CETS_CDP_NETWORK"), reason="set CETS_CDP_NETWORK=1 for the live portal test")
def test_live_portal(tmp_path):
    out = tmp_path / "live" / "10445.cets.json"
    r = _run(["to-cets", "portal:10445/TS_105_5@alignment=18924,voxel=4.99", "-o", str(out)])
    assert "reference_tomogram = 'TS_105_5_tomo_18956'" in r.stdout
    ds = load_dataset(out)
    recorded = run_data_from_json(json.loads((DATA / "10445_TS_105_5.json").read_text()))
    live_shift = ds.regions[0].alignments[0].projection_alignments[0].sequence[2].translation
    p = recorded.alignments[0].hub.per_section_alignment_parameters[0]
    assert live_shift == pytest.approx([p.x_offset * 1.54, p.y_offset * 1.54])


# ------------------------------------------------------------------ annotations


def test_recorded_annotations_to_cets(data, tmp_path):
    """Points: ``p = (loc - floor(N/2)) * s`` on the bound 4.99 A tomogram, matrices as stored; the mask keeps the
    tomogram grid read from the zarr metadata; the companion carries the ingestion metadata verbatim."""
    from cryoet_alignment.io.cets.annotations import annotation_kind, annotation_points

    ds, comp, sr = _to_cets(data, tmp_path / "10445.cets.json")
    region = ds.regions[0]
    kinds = {a.id: annotation_kind(a) for a in region.annotations}
    assert kinds == {
        "TS_105_5_tomo_18956_ann_69545_segmentationmask": "mask",
        "TS_105_5_tomo_18956_ann_69546_point": "points",
        "TS_105_5_tomo_18956_ann_69547_point": "points",
        "TS_105_5_tomo_18956_ann_175505_orientedpoint": "oriented_points",
    }
    for ann in region.annotations:
        if kinds[ann.id] == "mask":
            assert (ann.width, ann.height, ann.depth) == (1260, 1260, 368)
            assert ann.coordinate_transformations[0].sequence[1].scale == [4.99, 4.99, 4.99]
            assert ann.path.endswith("membrane-1.0_segmentationmask.zarr")
            assert comp.annotations[ann.id].mrc_path.endswith(".mrc") and comp.annotations[ann.id].mask_label == 1
            continue
        pid = int(ann.id.split("_ann_")[1].split("_")[0])
        nd = read_ndjson_points(
            DATA / "annotations" / f"{pid}_{'OrientedPoint' if kinds[ann.id] == 'oriented_points' else 'Point'}.ndjson"
        )
        r = annotation_points(ann, region)
        expect = (nd.locations - np.array([630.0, 630.0, 184.0])) * 4.99
        assert np.abs(r.points_a - expect).max() < 1e-9
        assert np.abs(r.points_corner_a / 4.99 - nd.locations).max() < 1e-9
        if nd.matrices is not None:
            assert np.abs(r.matrices - nd.matrices).max() < 1e-12  # the portal matrix IS the particle->tomogram matrix
        c = comp.annotations[ann.id]
        assert c.portal_annotation_id == pid and c.metadata["annotation_ingest_id"] and c.voxel_spacing_a == 4.99
    assert any(g.name == "mask_grid_matches_tomogram" and g.passed for g in sr.gates)
    assert "annotations" in [p["option"] for p in sr.provenance]


def test_annotations_staged_as_relion4_stars_and_config_blocks(data, tmp_path):
    out = tmp_path / "10445.cets.json"
    _to_cets(data, out)
    stage = tmp_path / "stage"
    r = _run(["from-cets", str(out), "-o", str(stage), "--deposition-id", "1", "--no-validate"])
    assert "[ok ] backend_point_parser_reproduces_points" in r.stdout
    assert "[ok ] backend_point_parser_reproduces_rotations" in r.stdout
    assert "mask file is not local" in r.stderr
    stars = sorted((stage / "annotations" / RUN).glob("*.star"))
    assert [p.name for p in stars] == [
        "10310_beta-amylase-1_point.star",
        "10310_ferritin-complex-1_point.star",
        "10358_apo-ferritin-octopi-1_orientedpoint.star",
    ]
    loc, mats = relion4_star_to_points(stars[2])  # the backend's own relion4 arithmetic
    nd = read_ndjson_points(DATA / "annotations" / "175505_OrientedPoint.ndjson")
    assert np.abs(loc - nd.locations).max() < 1e-6
    assert np.abs(np.einsum("nij,nkj->nik", mats, nd.matrices) - np.eye(3)).max() < 1e-6
    cfg = yaml.safe_load((stage / "ingestion_config.yaml").read_text())
    blocks = cfg["annotations"]
    assert len(blocks) == 4
    by_glob = {next(iter(b["sources"][0].values())).get("glob_string"): b for b in blocks}
    ori = by_glob["annotations/{run_name}/10358_apo-ferritin-octopi-1_orientedpoint.star"]
    assert (
        set(ori["sources"][0]) == {"OrientedPoint"}
        and ori["sources"][0]["OrientedPoint"]["file_format"] == "relion4_star"
    )
    assert ori["sources"][0]["OrientedPoint"]["order"] == "xyz" and "columns" not in ori["sources"][0]["OrientedPoint"]
    assert ori["metadata"]["annotation_ingest_id"] == "apo-ferritin-octopi-1" and "files" not in ori["metadata"]
    pt = by_glob["annotations/{run_name}/10310_ferritin-complex-1_point.star"]
    assert pt["sources"][0]["Point"]["columns"] == "xyz" and pt["metadata"]["annotation_object"]["id"] == "GO:0070288"
    mask = by_glob[TODO]
    assert (
        set(mask["sources"][0]) == {"SemanticSegmentationMask"}
        and mask["sources"][0]["SemanticSegmentationMask"]["mask_label"] == 1
    )
    assert "annotations: source glob is a TODO" in r.stdout + r.stderr
    # only one annotation, and none at all
    stage2 = tmp_path / "stage2"
    _run(
        [
            "from-cets",
            str(out),
            "-o",
            str(stage2),
            "--deposition-id",
            "1",
            "--no-validate",
            "--annotation",
            "TS_105_5_tomo_18956_ann_69547_point",
        ]
    )
    assert [p.name for p in (stage2 / "annotations" / RUN).glob("*")] == ["10310_beta-amylase-1_point.star"]
    stage3 = tmp_path / "stage3"
    _run(["from-cets", str(out), "-o", str(stage3), "--deposition-id", "1", "--no-validate", "--no-annotations"])
    assert not (stage3 / "annotations").exists() and "annotations" not in yaml.safe_load(
        (stage3 / "ingestion_config.yaml").read_text()
    )


@pytest.mark.skipif(not os.environ.get("CRYOET_DATA_PORTAL_BACKEND_PATH"), reason="set CRYOET_DATA_PORTAL_BACKEND_PATH")
def test_annotation_blocks_validate_against_backend_models(data, tmp_path):
    import sys

    out = tmp_path / "10445.cets.json"
    _to_cets(data, out)
    stage = tmp_path / "stage"
    _run(["from-cets", str(out), "-o", str(stage), "--deposition-id", "1", "--no-validate"])
    cfg = yaml.safe_load((stage / "ingestion_config.yaml").read_text())
    codegen = Path(os.environ["CRYOET_DATA_PORTAL_BACKEND_PATH"]) / "schema" / "ingestion_config" / "v1.0.0" / "codegen"
    sys.path.insert(0, str(codegen))
    try:
        from ingestion_config_models import AnnotationEntity  # type: ignore

        for block in cfg["annotations"]:
            AnnotationEntity.model_validate(block)
    finally:
        sys.path.remove(str(codegen))
