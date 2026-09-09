# cets-cryoet-data-portal

Convert cryoET Data Portal runs to the CETS cryo-ET standard, and stage CETS datasets for portal ingestion.

The package implements the rigid profile `cets-rigid/0.1` of CETS, documented in
[`cryoet-alignment/docs/cets.md`](https://github.com/uermel/cryoet-alignment/blob/uermel/cets/docs/cets.md).
The command is `cets-cdp`.

- `cets-cdp to-cets` reads runs through the portal API client and writes one CETS dataset. A whole
  portal dataset, single runs, or a specific alignment and reference tomogram can be requested.
- `cets-cdp from-cets` writes a staging directory in the shape portal ingestion configs expect, plus a
  draft `ingestion_config.yaml` with the derived values filled in and everything else marked for the
  curator. Nothing is inferred.

The portal's rigid per-section parameters (the `alignment_metadata.json` written at ingestion) are the
hub model of the codec, so portal to CETS to `.aln` reproduces the portal's own alignment file.

## Installation

```bash
pip install git+https://github.com/TomoBabel/cets-cryoet-data-portal.git
```

The package depends on `cryoet-alignment >= 0.3.0` (the CETS codec), `cryoet-data-portal >= 4.3.1` (the
API client) and the pinned `cets_data_model` commit listed in `pyproject.toml`. Until cryoet-alignment
0.3.0 is on PyPI, install it from its branch first:

```bash
pip install git+https://github.com/uermel/cryoet-alignment.git@uermel/cets
```

## Portal to CETS

### Source specification

`to-cets` takes portal references instead of files:

```
portal:<dataset>                              every run of the dataset
portal:<dataset>/<run>                        one run
portal:<dataset>/<run>@alignment=<id>         one run, a specific alignment
portal:<dataset>/<run>@alignment=<id>,voxel=<Å>  ... and the reference tomogram by voxel spacing
```

Without `@alignment=` the portal-standard alignment with the lowest id is used. Without `@voxel=` the
reference tomogram is the portal-standard one at the finest voxel spacing.

What is read per run:

```
TiltSeries          pixel_spacing, size, tilt_axis, acceleration_voltage, spherical_aberration_constant,
                    is_aligned (refused when true), https_mrc_file
Frames              acquisition_order, exposure_dose, accumulated_dose, https_frame_path
PerSectionParameters raw_angle, major/minor defocus, astigmatic_angle, phase_shift (radians -> degrees)
Alignment           alignment_metadata.json (per-section rotation, tilt, offsets in pixels), volume_dimension,
                    affine_transformation_matrix / volume_offset (refused unless identity / zero)
Tomogram            voxel_spacing, size, processing, reconstruction_method, ctf_corrected, https paths
FrameAcquisitionFile  the mdoc URL (companion, for the ingestion draft)
```

### Example

```bash
cets-cdp to-cets "portal:10445/TS_105_5@alignment=18924,voxel=4.99" -o cets/10445.cets.json
```

What happens per run:

1. Every raw section becomes a `TiltImage` with its nominal angle, pre-exposure dose and CTF; frames
   become movie stacks with `https://` (or `s3://`) paths.
2. Every portal tomogram becomes a `Tomogram` at the voxel spacing the portal declares.
3. The alignment box is checked against the reference tomogram's extent (gate
   `portal_volume_box_vs_reference_tomogram`) and the alignment is written as `portal<id>` in the
   centred physical frames, with AreTomo3's centre convention corrected on the way in.
4. Acquisition order, per-image exposure, kV/Cs, the tilt axis, `processing`, `reconstruction_method`
   and the mdoc URL go to the companion.

Console output (abridged):

```
== TS_105_5
   uri_scheme = 'https'  [default]
   pix = 1.54  [discovered]  (TiltSeries.pixel_spacing)
   reference_tomogram = 'TS_105_5_tomo_18956'  [discovered]  (portal tomogram 18956 at 4.99 Å (--voxel))
   voltage = 300.0  [discovered]  (TiltSeries.acceleration_voltage)
   cs = 2.7  [discovered]  (TiltSeries.spherical_aberration_constant)
   [ok ] portal_volume_box_vs_reference_tomogram value={'x': 6287.4, 'y': 6287.4, 'z': 1836.32} expected={...}
   [ok ] rows value=31
   WARNING: uri_scheme defaulted to 'https'; set it with --uri-scheme or config key series.TS_105_5.uri_scheme
   dropped: LOCAL alignment: the portal metadata carries the rigid per-section parameters only
wrote cets/10445.cets.json (1 region(s)) + 10445.cets-companion.json
report: cets/10445.cets.report.json
```

Output:

```
cets/
├── 10445.cets.json            # one Region per run: tilt series, movie stacks, alignment "portal<id>", tomograms
├── 10445.cets-companion.json  # acquisition order, exposure, kV/Cs, tilt axis, processing, method type,
│                              #   is_portal_standard, mdoc URL, header vs implied voxels, what was dropped
├── 10445.cets.report.json     # per run: provenance of every value, gates, warnings, errors
└── .portal/                   # cached alignment_metadata.json per run (--cache DIR to move it)
```

More examples:

```bash
# a whole dataset (one region per run; failing runs are reported, the rest continue)
cets-cdp to-cets portal:10445 -o cets/10445.cets.json

# several runs, s3 paths
cets-cdp to-cets portal:10445/TS_105_5 portal:10445/TS_106_1 -o cets/two.cets.json --uri-scheme s3

# the reference tomogram at a coarser voxel spacing than the portal-standard one
cets-cdp to-cets "portal:10445/TS_105_5@voxel=10.012" -o cets/10445_coarse.cets.json
```

## CETS to portal staging

```bash
cets-cdp from-cets cets/at3.cets.json -o staging/ --deposition-id 10301 \
    --method-type projection_matching --portal-standard --template curator.yaml
```

What happens per region:

1. The alignment and the reference tomogram are selected (flags required only when a region has
   several).
2. The alignment is staged in a format the backend parses: `.aln` when every volume X rotation is zero,
   otherwise IMOD `.xf/.tlt/.xtilt`. The staged files are read back with the same cryoet-alignment
   readers the backend's ingestion uses and must reproduce the document (gate
   `backend_parser_reproduces_hub`).
3. `.rawtlt` (nominal angles) and, when every image carries CTF metadata, a CTFFIND-parsable `_CTF.txt`
   are written.
4. Local tilt series, tomogram, mdoc and frame files named by the document or the companion are
   symlinked into the staging tree. Remote paths are not fetched.
5. The config draft is assembled: derived values (pixel spacing, tilt axis, tilt range and step,
   alignment block per format, CTF, rawtilt, voxel spacing) are filled; per-run values that differ go to
   `run_to_data_map.tsv`; everything else comes from `--template` or stays `TODO(curator)`.
6. Three checks run and are reported separately: every source glob resolves in the staging tree; the
   backend's generated schema (in-process); the backend's extended validator (in its conda env). The two
   backend checks are reported as blocked while placeholders remain, and skipped without the backend.

Console output (abridged):

```
== TS_105_5
   reference_tomogram = 'TS_105_5_tomo_18956'  [discovered]  (companion)
   voltage = 300.0  [companion]
   cs = 2.7  [companion]
   method_type = 'projection_matching'  [companion]
   portal_standard = True  [companion]
   [ok ] x_rotation value=0.0  staged as .aln
   [ok ] backend_parser_reproduces_hub value=1.1e-13 expected='< 1e-3 (px / deg; file precision)'
   [ok ] ctf_parses_with_one_header_line expected=31
   WARNING: no local tomogram file to stage: the tomograms block is left as a TODO
wrote staging/ingestion_config.yaml
TODO(curator): standardization_config.source_prefix
TODO(curator): datasets[0].metadata.dataset_title
...
validation: sources FAILED; schema skipped; extended skipped; 25 placeholder(s)
   tiltseries: source glob is a TODO
report: staging/cets_cdp.report.json
```

Output:

```
staging/
├── alignment/TS_01/TS_01.aln                 # every volume_x_rotation == 0 -> format ARETOMO3
├── alignment/TS_02/TS_02.{xf,tlt,xtilt}      # otherwise -> format IMOD
├── rawtlt/TS_01/TS_01.rawtlt                 # nominal_tilt_angle per raw section
├── ctf/TS_01/TS_01_CTF.txt                   # when every image carries CTF metadata
├── tiltseries/TS_01/TS_01.mrc -> ...         # symlink when the document's path is local
├── tomograms/TS_01/TS_01.mrc -> ...          # symlink when the tomogram path is local
├── collection_metadata/TS_01/TS_01.mdoc -> ...   # symlink when the companion names a local mdoc
├── frames/TS_01/... -> ...                   # symlinks when the movie stacks are local
├── run_to_data_map.tsv                       # per-run values when they differ between runs
├── ingestion_config.yaml                     # the draft
└── cets_cdp.report.json                      # provenance, gates, the three validation results
```

The draft, abridged. `float {pixel_spacing}` style references point into `run_to_data_map.tsv`:

```yaml
standardization_config: {deposition_id: 10301, source_prefix: TODO(curator)}   # + run_data_map_file when runs differ
alignments:
  - metadata: {format: ARETOMO3, alignment_type: GLOBAL, method_type: projection_matching, is_portal_standard: true}
    sources: [{source_glob: {list_glob: 'alignment/{run_name}/*.aln'}}]
tiltseries:
  - metadata: {pixel_spacing: 1.54, tilt_axis: -96.33, tilt_range: {min: -45.0, max: 45.0}, tilt_step: 3.0,
               is_aligned: false, acceleration_voltage: 300000, camera: TODO(curator), ...}
    sources: [{source_glob: {list_glob: TODO(curator)}}]
rawtilts:
  - sources: [{source_glob: {list_glob: 'rawtlt/{run_name}/*.rawtlt'}}]
ctfs:
  - metadata: {format: CTFFIND}
    sources: [{source_glob: {list_glob: 'ctf/{run_name}/*_CTF.txt'}}]
```

The backend rules the draft follows: a `tiltseries` block needs a `collection_metadata` block, which in
turn needs a `frames` block. When no mdoc is known the `collection_metadata` block is left out with a
warning, and when no frames are staged the `frames` block is a literal `default`.

More examples:

```bash
# a Warp project converted by cets-warpm, tilt series already in the bucket
cets-cdp from-cets cets/warp.cets.json -o staging/ --deposition-id 10302 \
    --tiltseries-glob 'TiltSeries/{run_name}/*.mrc' --template curator.yaml

# skip the backend validators (e.g. on a machine without the backend checkout)
cets-cdp from-cets cets/at3.cets.json -o staging/ --deposition-id 10301 --no-validate

# run the validators
export CRYOET_DATA_PORTAL_BACKEND_PATH=/path/to/cryoet-data-portal-backend
export CRYOET_DATA_PORTAL_BACKEND_CONDA_ENV=ingestion
cets-cdp from-cets cets/at3.cets.json -o staging/ --deposition-id 10301 --template curator.yaml
```

A template is an ordinary ingestion config written by the curator. Its `standardization_config.source_prefix`,
`datasets`, `tiltseries` metadata, `tomograms` metadata and `frames` metadata are merged into the draft;
derived values always win.

## Command reference

### `cets-cdp to-cets`

```
cets-cdp to-cets [OPTIONS] SOURCES...
```

`SOURCES` are portal references (above).

| Option | Meaning | Derived from, when absent |
|---|---|---|
| `-o, --output FILE` | CETS dataset JSON to write (required) | |
| `--name TEXT` | dataset name | the portal dataset id |
| `--uri-scheme https\|s3` | scheme of the paths written into the document | `https`, with a warning |
| `--cache DIR` | cache for fetched `alignment_metadata.json` | `OUT_DIR/.portal` |
| `--voltage kV`, `--cs mm`, `--amp-contrast F` | companion values only | portal kV / Cs; amplitude contrast absent |
| `--fail-fast` | stop at the first failing run | continue, exit 1 at the end |
| `--overwrite` | replace existing outputs | error when outputs exist |
| `--config FILE` | YAML config with overrides (see below) | |

### `cets-cdp from-cets`

```
cets-cdp from-cets [OPTIONS] DOCUMENT
```

`DOCUMENT` is a CETS dataset JSON. Its companion manifest is read when present.

| Option | Meaning | Derived from, when absent |
|---|---|---|
| `-o, --output DIR` | staging directory (required) | |
| `--deposition-id N` | the portal deposition (required) | |
| `--region ID` | region(s) to convert; repeatable | all |
| `--alignment NAME\|N` | alignment instance name or 0-based index | required when a region has several |
| `--tomogram ID` | reference tomogram id | companion; required when a region has several |
| `--method-type` | `fiducial_based`, `patch_tracking`, `projection_matching`, `simulated`, `undefined` | companion; else `TODO(curator)` |
| `--portal-standard / --no-portal-standard` | `is_portal_standard` of the alignment block | companion; else `TODO(curator)` |
| `--template FILE` | curator-written config to merge (see above) | placeholders |
| `--tiltseries-glob PATTERN` | config glob for tilt series that are not staged | symlink + glob when local; else `TODO(curator)` |
| `--tomograms-glob PATTERN` | config glob for tomograms that are not staged | symlink + glob when local; else the block is a `TODO` |
| `--voltage kV`, `--cs mm` | tiltseries metadata | companion |
| `--defocus-hand -1\|1` | defocus handedness of the staged `_CTF.txt` | companion; column omitted otherwise |
| `--no-ctf` | do not stage CTFs | off |
| `--validate / --no-validate` | run the backend validators when the env vars are set | on |
| `--fail-fast`, `--overwrite`, `--config FILE` | as above | |

## Values, defaults and the config file

Every value a command needs resolves through one chain:

```
CLI flag  >  --config FILE  >  companion manifest  >  discovered (portal API)  >  package default
```

Each resolution is printed as `option = value  [source]` and recorded in the report. Falling back to a
package default always prints a warning that names the flag and the config key.

```yaml
cets:                          # every package and command
  voltage: 300
  cs: 2.7
cets-cdp:
  to-cets: {uri_scheme: s3}
  from-cets: {method_type: projection_matching, portal_standard: true, template: curator.yaml}
series:                        # per run; wins over the command section
  TS_01: {tiltseries_glob: 'TiltSeries/TS_01/*.mrc'}
```

```bash
cets-cdp from-cets cets/at3.cets.json -o staging/ --deposition-id 10301 --config cets.yaml
```

## What is refused, what is dropped

- Tilt series flagged `is_aligned`: the per-section shifts do not apply to an aligned stack.
- Alignments with a non-identity `affine_transformation_matrix` or a non-zero `volume_offset`.
- LOCAL alignments contribute their rigid per-section parameters; the local part is noted as dropped.
- Nothing in the ingestion draft is inferred: `processing`, `reconstruction_method`, camera, microscope,
  software and the dataset metadata come from the companion or the template, or stay placeholders.

## Development

```bash
pip install -e '.[dev]'
pytest                                   # recorded metadata of 10445/TS_105_5
CETS_CDP_NETWORK=1 pytest                # plus the live portal test
CRYOET_DATA_PORTAL_BACKEND_PATH=... pytest    # plus the schema check against the backend
pre-commit run --all-files               # ruff 0.11.12, ruff-format, mypy 1.8.0
```
