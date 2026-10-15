# CLAUDE.md

Guidance for Claude Code (claude.ai/code) in this repository. People start at the
[README](README.md); how the pipeline works is in [docs/architecture.md](docs/architecture.md).

## Project

FOMO is a Nextflow DSL2 pipeline that annotates candidate lncRNAs on target assemblies by
transferring the lncRNA annotation of source species: preprocessing per source, one alignment of
all sources per target, benchmarking, curation of the intron chains several species share, one
MultiQC report per target and one per run. Input: a samplesheet (`species,role,fasta,gff3`).
Output: [docs/output.md](docs/output.md).

## Developing

1. **Look for prior art first.** Search nf-core modules with the Seqera MCP tools; then vendor
   the module into `modules/local/`, trimmed to what FOMO calls (there is no nf-core module
   directory and no `modules.json`).
2. Keep the conventions below: meta maps, `ext.args`/`ext.prefix` in `conf/modules.config`, one
   `img_*` label per module.
3. A new image is a Wave (or BioContainers) image, pinned, declared once in
   `conf/containers.config`.
4. Before launching: `nextflow config -profile <p>`, `dev/nextflow_lint.sh`, the stub matrix.
5. Design decisions go in `plans/` (see [plans/README.md](plans/README.md)); history goes in the
   CHANGELOG and commit messages, never in code comments.

## Glossary

- **source / target / both** — the samplesheet roles: [README](README.md#samplesheet).
- **track** — a feature type plus the decoy flag: `lnc_RNA`, `lnc_RNA.decoy`, `mRNA`,
  `mRNA.decoy`. mRNA (`--include_mrna`) is a positive control, decoys (`--include_decoy`) a
  false-positive baseline; lncRNA is the deliverable.
- **allModels** — every source's projected models of a track, pooled (`allModels_raw` as a model
  id in reports and gffcompare files).
- **record** — a statistics record (`fomo-data/1`): [docs/data_records.md](docs/data_records.md).
- **S, T, F, D** — sources, targets, enabled feature types, 2 with decoys else 1.

## Invariants

Each line: the rule, then what enforces it.

1. The filtered annotation lists exactly the projected FASTA's records, matched by exact ID
   through the rename map; any disagreement fails. → `tests/filter_gff_by_id`, golden runs.
2. A record's species travels from `>tid|type|species` to the BAM read name, to `source=`, to the
   split; projected transcript IDs are `<tid>|<source>`, and a duplicate ID fails. →
   `tests/bam_to_gff`, `tests/gff_by_source`.
3. A species is never aligned to its own genome, on any track. → `tests/fasta_drop_species`, the
   golden runs (no `from_<target>` file).
4. Genomes and annotations stay compressed; there is no gunzip task and no `samtools faidx`.
   → golden V4 (plain gzip) must equal V1.
5. Decoys come from the raw annotation's gene-free space and are reproducible with the seed. →
   `tests/relocate_loci` (byte-stable fixture), golden V2.
6. The enabled feature types are derived once (`workflows/fomo.nf`); projection ⋈ reference
   pairs on `[target, feature type]`. → stub matrix task counts.
7. A target onto which nothing projects is reported with zero models. → `tests/gff_by_source`
   (empty manifest), the `remainder: true` join in `projection.nf`.
8. TD2 with an empty FASTA or no complete ORF means "nothing filtered". →
   `tests/modules/test_td2_predict.py`.
9. Transcript selection (stranded, spliced, longest isoform, ties to the smallest ID) happens only
   in `bin/filter_transcript.py`. → `tests/filter_transcript`.
10. Curated IDs are content hashes, stable across reruns, input order and source sets. →
    `tests/curate_models`.
11. One MultiQC report per target plus one aggregate-only run summary; the publish layout. →
    stub matrix (`tests/pipeline/expected/*.paths`), `tests/test_multiqc_registry.py`.
12. Task counts follow the formulas (`F·S` preprocessing, `F·D·T` alignments, ...). → stub
    matrix.
13. Reports read only records, and records carry every reported value; the run summary is
    byte-identical on the golden fixtures. → `tests/records`, `tests/run_summary_tables`.
14. Every module: one `img_*` label, a stub, versions on the topic (not MULTIQC), the tag rule,
    no `when:`; every `bin/` script with a shebang is executable; the vendored
    `bin/gff-feature-stats` matches its checksum. → `tests/test_repo_conventions.py`.
15. A process-level `ext.*` in `conf/modules.config` stays after every `withName` block (the
    config parser would otherwise capture the later blocks' `ext.*`). →
    `tests/test_repo_conventions.py`.
16. Deliverables never change in a refactor: published GFF3/BAM/gffcompare/curated files and
    `summary/tables`. → golden gate `tests/golden/gate.sh` against `tests/golden/V*.tsv`.

## Architecture map

Subworkflows (`subworkflows/local/`), wired by `workflows/fomo.nf`, never nested:

| subworkflow | emits |
|:---|:---|
| PREPROCESSING | `spliced_fasta`, `decoy_spliced_fasta`, `filtered_gff3`, `data` |
| PROJECTION | `allmodels_raw`, `gff3` (batch + manifest), `mqc` (samtools stats), `data` |
| BENCHMARKING | `mqc` (gffcompare stats), `data` |
| CURATION | `data` |
| REPORTING | — (one MULTIQC per target, rendering records) |
| RUN_SUMMARY | — (RUN_SUMMARY_TABLES + one MULTIQC) |

Meta keys: `id`, `role`, `target_id` (routes reports, names publish dirs), `feature_type`,
`decoy`, `track`/`label`/`order`/`td2` (stamped by `withTrack`), `kind`, `stage`,
`exclude_self`. Helpers: `subworkflows/local/fomo_utils.nf` (`withTrack`, `pairBy`,
`routeByTarget`, `softwareVersionsYaml`). Records: `bin/fomo_stats.py`; report sections:
`bin/fomo_sections.py`; shared I/O: `bin/fomo_io.py`. Publish layout: `docs/output.md`.

## Conventions

- **Comments** say what is true now and why, in at most ~3 lines. No run names, job IDs, clades,
  measurements, incident stories or history; `Rationale: plans/NN` only when that plan is the
  current rationale, and never in user-facing text. `dev/code_metrics.py --check` enforces it.
- **Naming.** Processes `VERB_OBJECT`, module file = process name in lower case; emits `mqc`
  (MultiQC tool files) and `data` (records); Sequence Ontology names (`lnc_RNA`, `mRNA`) in
  identifiers and file names, "lncRNA" in prose and labels.
- **Modules** carry one `img_*` label and a resource label, a stub, and the tag
  `"${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"`; no `container`,
  `conda` or `when:` lines. Optional branches are gated by an empty input channel.
- **Config.** `nextflow.config` includes `base.config` (every resource), `containers.config`,
  the profiles, then `modules.config` (`ext.*`, `publishDir`; selectors `.*:NAME$`). Site
  profiles (`conf/crg.config`) hold site settings only. Every param is also in
  `nextflow_schema.json` (unknown params abort the run).
- **Python** (`bin/`): stdlib only; `fomo_io` for files, errors (`FomoError` → exit 1) and
  logging; inputs `--gff/--fasta/--json`, outputs `--output` or `--prefix`; complexity ≤ 15;
  `ruff check` and `ruff format`.
- **Shell**: `set -euo pipefail`, `export LC_ALL=C`, shellcheck-clean; awk without `var (expr)`
  concatenation (busybox awk reads it as a call).
- **Containers.** Tools without Python write records with bash+awk (`bin/*_record.sh`); Python
  scripts run on `img_python`, `img_curate` or the MultiQC image.

## Checklists

**Add a track** (a new feature type): extend `trackOf` in `fomo_utils.nf`; derive it in
`workflows/fomo.nf`; check the `--feature-type` paths of `filter_transcript.py` and
`relocate_loci.py`; add rows to the stub matrix counts and a golden variant if it is a
deliverable.

**Add a statistics section:** write a record (a `SECTIONS` entry in `fomo_stats.py`, the
producer's `data` emit, a `docs/data_records.md` entry); render it (a `Section` in
`fomo_sections.py`, a `custom_data` block and `sp` pattern in `assets/multiqc/sections.yml`, a
place in `main.yml`'s `report_section_order`); test it in `tests/records`.

**Add a parameter:** `nextflow.config` `params`, `nextflow_schema.json`, the module's `ext.args`
in `conf/modules.config` (modules never read `params`), `docs/usage.md`.

**Add a module:** `modules/local/<name>.nf` with one `img_*` label, a resource label, the tag
rule, versions on the topic and a stub; its image in `conf/containers.config`; `ext.*` and
`publishDir` in `conf/modules.config`; pytest for any new `bin/` script; update the stub matrix
counts.

## Testing

| gate | command |
|:---|:---|
| unit | `pytest -q` (environment: `dev/environment.yml`) |
| lint | `pre-commit run --all-files` (ruff, ruff format, shellcheck, `dev/nextflow_lint.sh`, `dev/code_metrics.py --check`, `dev/check_links.py`) |
| stub matrix | `nf-test test tests/pipeline --profile +singularity` |
| golden | `tests/golden/gate.sh <scratch dir>` (V1–V4 on SLURM, compared with `tests/golden/`) |

CI (`.github/workflows/ci.yml`) runs lint, pytest, the stub matrix and the golden variants with
Docker. Never regenerate a golden file to make a gate pass.
