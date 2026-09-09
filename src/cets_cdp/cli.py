"""``cets-cdp`` command line: ``to-cets`` (portal API -> CETS) and ``from-cets`` (CETS -> staging + config draft)."""

import os
import sys
from pathlib import Path

import click
from cryoet_alignment.io.cets.cli_support import (
    Report,
    SeriesReport,
    common_options,
    echo,
    finish,
    load_config,
    make_resolver,
    parse_alignment_selector,
    print_series,
    selection_options,
    warn,
)
from cryoet_alignment.io.cets.companion import Companion
from cryoet_alignment.io.cets.entities import dataset_entity, dump_json, load_dataset, validate_document

from cets_cdp import __version__
from cets_cdp.api import fetch_run, find_runs, parse_portal_source, portal_client
from cets_cdp.from_cets import (
    build_config,
    check_extended,
    check_schema,
    check_sources_resolve,
    stage_run,
    todos,
    write_config,
)
from cets_cdp.to_cets import portal_to_cets

PACKAGE = "cets-cdp"
TO_CETS_OPTIONS = {"uri_scheme", "voltage", "cs", "amp_contrast"}
FROM_CETS_OPTIONS = {"method_type", "portal_standard", "template", "voltage", "cs", "defocus_hand", "no_ctf", "validate", "tiltseries_glob", "tomograms_glob"}


@click.group()
@click.version_option(__version__)
def main():
    """cryoET Data Portal <-> CETS (rigid profile cets-rigid/0.1)."""


@main.command("to-cets")
@click.argument("sources", nargs=-1, required=True)
@click.option("-o", "--output", "output", required=True, type=click.Path(dir_okay=False), help="CETS dataset JSON to write.")
@click.option("--name", default=None, help="Dataset name (default: the dataset id).")
@click.option("--uri-scheme", "uri_scheme", type=click.Choice(["https", "s3"]), default=None)
@click.option("--cache", "cache", type=click.Path(file_okay=False), default=None, help="Cache for fetched metadata [OUT_DIR/.portal].")
@click.option("--voltage", type=float, default=None, help="kV (companion only; default: the portal value).")
@click.option("--cs", type=float, default=None)
@click.option("--amp-contrast", "amp_contrast", type=float, default=None)
@common_options
def to_cets(sources, output, name, uri_scheme, cache, config_path, overwrite, fail_fast, **cli):
    """Convert portal runs (portal:<dataset>[/<run>][@alignment=ID,voxel=A]) to a CETS dataset JSON."""
    out = Path(output)
    if out.exists() and not overwrite:
        raise click.ClickException(f"{out} exists (use --overwrite)")
    out.parent.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path, TO_CETS_OPTIONS)
    flags = {k: v for k, v in cli.items() if v is not None}
    if uri_scheme:
        flags["uri_scheme"] = uri_scheme
    cache_dir = Path(cache) if cache else out.parent / ".portal"
    try:
        specs = [parse_portal_source(t) for t in sources]
    except ValueError as e:
        raise click.ClickException(str(e)) from e
    client = portal_client()
    report = Report(PACKAGE, "to-cets")
    companion = Companion(generator=f"{PACKAGE} {__version__}")
    regions = []
    dataset_ids = set()
    for spec in specs:
        for run in find_runs(client, spec):
            sr = SeriesReport(run.name)
            report.series.append(sr)
            try:
                data = fetch_run(client, run, cache_dir=cache_dir)
                res = make_resolver(PACKAGE, "to-cets", flags, config, run.name, sr)
                result = portal_to_cets(data, res, sr, alignment_id=spec.alignment_id, voxel=spec.voxel_spacing)
            except Exception as e:  # noqa: BLE001
                sr.error = str(e)
                print_series(sr)
                if fail_fast:
                    break
                continue
            dataset_ids.add(data.dataset_id)
            regions.append(result.region)
            companion.tilt_series[run.name] = result.tilt_series_companion
            if result.alignment_companion is not None:
                companion.alignments.append(result.alignment_companion)
            companion.tomograms.update(result.tomogram_companions)
            print_series(sr)
    if regions:
        ds = dataset_entity(name or ("-".join(str(d) for d in sorted(dataset_ids))), regions)
        validate_document(ds)
        dump_json(ds, out)
        companion.dump(Companion.path_for(out))
        echo(f"wrote {out} ({len(regions)} region(s)) + {Companion.path_for(out).name}")
    finish(report, out.with_name(out.name.split(".")[0] + ".cets.report.json"))


@main.command("from-cets")
@click.argument("document", type=click.Path(exists=True, dir_okay=False))
@click.option("-o", "--output", "output", required=True, type=click.Path(file_okay=False), help="Staging directory.")
@click.option("--deposition-id", "deposition_id", type=int, required=True)
@selection_options
@click.option("--method-type", "method_type", type=click.Choice(["fiducial_based", "patch_tracking", "projection_matching", "simulated", "undefined"]), default=None)
@click.option("--portal-standard/--no-portal-standard", "portal_standard", default=None)
@click.option("--template", type=click.Path(exists=True, dir_okay=False), default=None, help="Curator-written config whose dataset/deposition/tiltseries/tomograms metadata is merged in.")
@click.option("--tiltseries-glob", "tiltseries_glob", default=None, help="Config glob for tilt series that are not staged (e.g. already in the bucket).")
@click.option("--tomograms-glob", "tomograms_glob", default=None, help="Config glob for tomograms that are not staged.")
@click.option("--voltage", type=float, default=None)
@click.option("--cs", type=float, default=None)
@click.option("--defocus-hand", "defocus_hand", type=click.Choice(["-1", "1"]), default=None)
@click.option("--no-ctf", "no_ctf", is_flag=True, default=None)
@click.option("--validate/--no-validate", "validate", default=None, help="Run the backend validators when the backend env vars are set [on].")
@common_options
def from_cets(document, output, deposition_id, regions, alignment, tomogram, config_path, overwrite, fail_fast, **cli):
    """Stage alignments/rawtilts/CTFs (+ local data links) and write an ingestion-config draft."""
    import yaml

    doc = Path(document)
    ds = load_dataset(doc)
    companion = Companion.load_for(doc)
    staging = Path(output)
    staging.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path, FROM_CETS_OPTIONS)
    flags = {k: v for k, v in cli.items() if v is not None}
    if "defocus_hand" in flags:
        flags["defocus_hand"] = int(flags["defocus_hand"])
    report = Report(PACKAGE, "from-cets")
    wanted = set(regions)
    per_run = []
    for region in ds.regions:
        if wanted and region.id not in wanted:
            continue
        sr = SeriesReport(region.id)
        report.series.append(sr)
        try:
            res = make_resolver(PACKAGE, "from-cets", flags, config, region.id, sr)
            per_run.append(stage_run(
                region, res, sr, staging=staging, doc_dir=doc.parent, companion=companion,
                alignment_selector=parse_alignment_selector(alignment), tomogram_selector=tomogram, overwrite=overwrite,
            ))
        except Exception as e:  # noqa: BLE001
            sr.error = str(e)
        print_series(sr)
        if sr.error and fail_fast:
            break
    extra = {}
    if per_run:
        template = None
        tpl_path = flags.get("template")
        if tpl_path is None and config is not None:
            hit = config.lookup("template", PACKAGE, "from-cets")
            tpl_path = hit if isinstance(hit, str) else None
        if tpl_path:
            template = yaml.safe_load(Path(tpl_path).read_text())
        cfg, rows = build_config(per_run, deposition_id=deposition_id, template=template, staging=staging)
        cfg_path = write_config(cfg, rows, staging)
        echo(f"wrote {cfg_path}" + (" + run_to_data_map.tsv" if rows else ""))
        todo = todos(cfg)
        for t in todo:
            warn(f"TODO(curator): {t}")
        problems = check_sources_resolve(cfg, staging, per_run, rows)
        do_validate = flags.get("validate", True)
        backend = os.environ.get("CRYOET_DATA_PORTAL_BACKEND_PATH")
        env = os.environ.get("CRYOET_DATA_PORTAL_BACKEND_CONDA_ENV")
        if not do_validate:
            schema_status, schema_msgs = "skipped", ["--no-validate"]
            ext_status, ext_msgs = "skipped", ["--no-validate"]
        elif todo:
            schema_status, schema_msgs = "blocked", [f"{len(todo)} placeholder(s) remain; fill them (or pass --template) first"]
            ext_status, ext_msgs = "blocked", schema_msgs
        else:
            schema_status, schema_msgs = check_schema(cfg, backend)
            ext_status, ext_msgs = check_extended(cfg_path, backend, env, staging / "validation")
        extra["validation"] = {
            "sources_resolve": {"status": "passed" if not problems else "failed", "problems": problems},
            "schema": {"status": schema_status, "messages": schema_msgs},
            "extended": {"status": ext_status, "messages": ext_msgs},
            "placeholders": todo,
        }
        echo(f"validation: sources {'ok' if not problems else 'FAILED'}; schema {schema_status}; extended {ext_status}; {len(todo)} placeholder(s)")
        for p in problems:
            warn(f"   {p}")
        if schema_status == "failed":
            warn("   schema: " + schema_msgs[0][:500])
        report.extra = extra
        real_problems = [p for p in problems if "is a TODO" not in p]
        if real_problems or schema_status == "failed" or ext_status == "failed":
            from cryoet_alignment.io.cets.cli_support import Gate

            report.series[0].gates.append(Gate("config_validation", False, note="see report 'validation'"))
    finish(report, staging / "cets_cdp.report.json")


if __name__ == "__main__":
    sys.exit(main())
