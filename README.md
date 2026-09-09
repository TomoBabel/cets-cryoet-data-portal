# cets-cryoet-data-portal

Converters between the cryoET Data Portal and the CETS cryo-ET standard, rigid profile `cets-rigid/0.1`
(see `cryoet-alignment/docs/cets.md`). Command name: `cets-cdp`.

```
cets-cdp to-cets portal:10445/TS_105_5@alignment=18924,voxel=4.99 -o out/10445.cets.json
cets-cdp from-cets out/run.cets.json -o staging/ --deposition-id 10301 --method-type projection_matching \
    --portal-standard --template curator.yaml
```

`to-cets` reads a run through the API client: tilt series (pixel spacing, size, tilt axis, kV, Cs,
`is_aligned` — refused when true), frames (acquisition order, exposure, accumulated dose, frame URLs),
per-section parameters (nominal angles, CTF; phase shift in radians converted to degrees), the alignment's
`alignment_metadata.json` (the same per-section parameters the ingestion wrote; offsets are pixels) and the
tomograms. Tomograms carry the voxel spacing the portal declares, and the reference frame of the alignment
is the portal-standard tomogram at the finest voxel (or the one selected with `@voxel=`), so the alignment
box equals the portal's own `volume_dimension`. The companion records, for information only, the voxel the
raw field of view would imply (`pixel_spacing × size_x_ts / size_x_tomo`). Non-identity registration
matrices / offsets are refused; LOCAL alignments contribute their rigid part with a note. Paths are
`https://` (or `s3://` with `--uri-scheme s3`).

`from-cets` stages, per region, `alignment/<run>/<run>.aln` (or `.xf/.tlt/.xtilt` when the alignment
carries X rotations — the backend reads that combination), `rawtlt/<run>/<run>.rawtlt`,
`ctf/<run>/<run>_CTF.txt` (CTFFIND-parsable), links for local tilt series / tomogram / mdoc / frame files,
`run_to_data_map.tsv` when runs differ, and `ingestion_config.yaml`: derived values filled (pixel spacing,
tilt axis = median rotation, tilt range/step, alignment block(s), CTF block, voxel spacing), everything
else from `--template` or left as `TODO(curator)` — nothing is inferred. Three checks are reported
separately: every source glob resolves in the staging directory; the backend's generated schema
(in-process, `CRYOET_DATA_PORTAL_BACKEND_PATH`); the backend's extended validator in its conda env
(`CRYOET_DATA_PORTAL_BACKEND_CONDA_ENV`). Schema checks are blocked while placeholders remain.

| option | default / derivation |
|---|---|
| `--uri-scheme https\|s3` | `https` (warned) |
| `--voltage --cs --amp-contrast` | portal kV / Cs; companion only |
| from-cets `--method-type`, `--portal-standard` | companion (portal `alignment_method`, `is_portal_standard`) else `TODO(curator)` |
| from-cets `--template config.yaml` | curator blocks merged (datasets, tiltseries, tomograms, frames metadata) |
| from-cets `--tiltseries-glob`, `--tomograms-glob` | config globs for data that is not staged |
| from-cets `--alignment NAME\|N`, `--tomogram ID` | required when a region has several |
| from-cets `--validate/--no-validate` | on |
