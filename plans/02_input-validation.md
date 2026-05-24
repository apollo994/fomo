# Input validation with nf-schema

Replace the hand-written samplesheet parser in `workflows/fomo.nf` with the
nf-core standard `nf-schema` plugin.

Status: **implemented** on 2026-05-24. End-to-end `test,docker` run passes
(128/128 tasks).

## Goal

- Declarative samplesheet validation (required fields, enum on `role`, file
  existence, extension patterns).
- Single source of truth for samplesheet structure (`assets/schema_input.json`).
- Auto-built meta map; no more manual `assert` + path resolution.
- Foundation for future `validateParameters()` on pipeline params.

## Scope

Samplesheet validation only. Params-level validation
(`nextflow_schema.json` + `validateParameters()`) is out of scope here — can
follow in a separate plan once params stabilise.

## Files changed

1. **`nextflow.config`** — plugin declaration.
   ```groovy
   plugins {
       id 'nf-schema@2.3.0'
   }
   ```

2. **`assets/schema_input.json`** — new. Uses JSON-schema **draft 2020-12**
   (required by nf-schema 2.x; draft-07 is rejected). Top-level shape is an
   `array` of `items`, as required by `samplesheetToList`. Columns:
   - `species` → string, required, pattern `^\S+$`, `meta: ["id"]`
   - `role` → string, required, `enum: ["source", "target"]`, `meta: ["role"]`
   - `fasta` → string, required, `format: "file-path"`, `exists: true`,
     pattern `^\S+\.(fa|fasta|fna)(\.gz)?$`
   - `gff3` → string, required, `format: "file-path"`, `exists: true`,
     pattern `^\S+\.gff3?(\.gz)?$`

   Each property also has an `errorMessage` for nicer reporting.

3. **`workflows/fomo.nf`** — manual parsing block (22 lines) replaced with:
   ```groovy
   include { samplesheetToList } from 'plugin/nf-schema'

   Channel
       .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
       .branch { meta, fasta, gff3 ->
           target: meta.role == 'target'
           source: meta.role == 'source'
       }
       .set { ch_input }
   ```
   Downstream `PREPROCESSING(ch_input.source)` and `STATS(...)` calls are
   unchanged — tuple shape `(meta, fasta, gff3)` is preserved.

4. **`assets/samplesheet.csv`** — **unchanged**. Paths stay relative to the
   launch directory (e.g. `assets/test_data/...`). See "Path resolution" below.

## Path resolution — important correction

The initial plan claimed nf-schema resolves paths relative to the samplesheet
file. **This is wrong.** `samplesheetToList` resolves relative paths against
the **launch directory (CWD)** — the same place Nextflow itself resolves
relative paths. Absolute paths and URIs (S3, etc.) are passed through
unchanged.

Implication: as long as the user runs `nextflow run .` from the repo root,
paths like `assets/test_data/...` work. For portability across launch
directories, future samplesheets should prefer absolute paths or remote URIs
(the standard nf-core convention).

## Behaviour changes vs. the old parser

- **`projectDir` fallback removed**: the old code prefixed `${projectDir}/` to
  any non-absolute path. nf-schema does not do this — relative paths now
  resolve against CWD. For the bundled `test` profile this is equivalent in
  practice because the user launches from the repo root, but downstream
  users with samplesheets outside the repo will see different behaviour.
- **Error reporting**: nf-schema aggregates all row errors before failing,
  instead of stopping on the first `assert`. Verified during testing: a
  samplesheet with bad paths reported all 10 errors (5 rows × 2 fields) in
  one shot.

## Validation evidence

- `nextflow run . -profile test,docker` → 128/128 tasks succeeded, ~1 min.
- Earlier failing runs proved validation works: a draft-07 schema was
  rejected with a clear migration pointer; a samplesheet with wrong relative
  paths produced one error line per (row, field).

## Issues encountered during implementation

- **draft-07 → draft 2020-12**: nf-schema 2.x mandates 2020-12. First run
  failed with a clear error and migration link.
- **Misunderstood path resolution**: assumed samplesheet-relative; actually
  CWD-relative. Reverted the samplesheet path edit after the first failing
  run made this obvious.

## Out of scope (follow-ups)

- `nextflow_schema.json` + `validateParameters()` for `params.input`,
  `params.outdir`, `params.relocate_seed`, etc.
- `paramsSummaryLog()` / `paramsSummaryMap()` to print resolved params at
  workflow start.
- Switching the bundled samplesheet to absolute paths / remote URIs so the
  pipeline runs cleanly from any launch directory.
