# Plan: refactor STATS → REPORTING, stage-local QC, split MultiQC config

## Context

Today `subworkflows/local/stats.nf` does two things in one place:

1. **Stat production** — runs AGAT on raw/filtered/decoy GFFs and SeqKit on
   source/decoy FASTAs, then runs custom adapter modules (`AGAT_TO_MQC`,
   `SEQKIT_TO_MQC`) to emit MultiQC-shaped TSVs.
2. **Report aggregation** — calls MULTIQC over the union of everything.

This works for the preprocessing + benchmarking surface we have today, but:

- `STATS.take:` is up to **7 channels** — every new stat producer adds another
  channel and another piece of upstream knowledge to STATS.
- `assets/multiqc_config.yml` is monolithic and `extra_fn_clean_exts` is a
  hidden coupling between every filename in the pipeline and report.
- Future stages (alignment QC, projection QC, splice-junction validation)
  will each add more channels and more config sections.

This plan refactors the boundary so that **stat production lives in the
subworkflow that produces the artifacts**, REPORTING becomes a thin
assembler, and the MultiQC config splits per stage. Functionality is
unchanged — only the seams move.

## Decisions

| Decision | Choice |
|----------|--------|
| Subworkflow naming | `STATS` → `REPORTING`. The new name describes what it does (run MultiQC), not what flows through it. |
| Where AGAT/SeqKit run | Inside `PREPROCESSING` (it owns source raw + filtered + decoy + spliced FASTAs). |
| Where target raw GFF stats run | Inside `BENCHMARKING` — it already owns target-GFF decompression (GUNZIP_TARGET_GFF) and filtering (FILTER_TARGET). |
| `mqc_files` channel type | `Channel<tuple(meta, path)>` — same shape as the rest of the pipeline. REPORTING strips meta before handing to MULTIQC. Keeps traceability (e.g. logging which subworkflow emitted which file) at the cost of one `.map` per emit. |
| Empty `mqc_files` | Subworkflows with no current stats (`PROJECTION`) emit `Channel.empty()` so the top-level mix is uniform. |
| Multiple MultiQC configs | Split into **two** YAML files and use MULTIQC's existing two config slots (`multiqc_config` + `extra_multiqc_config`). No new merge module needed. |
| Report shape | Stays single report (`multiqc_report.html`). Per-stage reports are deferred until the single report gets too tall — the config split is the prerequisite. |
| Module-config selector renames | `FOMO:STATS:*` → `FOMO:REPORTING:*` (just MULTIQC remains). Stat-producer selectors move to `FOMO:PREPROCESSING:*` / `FOMO:BENCHMARKING:*`. |

## 1. New `REPORTING` subworkflow

Replace `subworkflows/local/stats.nf` with `subworkflows/local/reporting.nf`:

```groovy
include { MULTIQC } from '../../modules/nf-core/multiqc/main'

workflow REPORTING {
    take:
    ch_mqc_files   // Channel<tuple(meta, path)> — *_mqc.tsv / *.stats files

    main:

    ch_files = ch_mqc_files.map { _meta, p -> p }.collect()

    ch_mqc_input = ch_files.map { files ->
        tuple(
            [id: 'multiqc'],
            files,
            file(params.multiqc_main_config,     checkIfExists: true),
            file(params.multiqc_sections_config, checkIfExists: true),
            [],
            []
        )
    }

    MULTIQC(ch_mqc_input)

    emit:
    report = MULTIQC.out.report
    data   = MULTIQC.out.data
}
```

No new local module needed — the existing nf-core MULTIQC signature takes
two config slots (`multiqc_config` at position 3, `extra_multiqc_config` at
position 4), and we fill both.

## 3. Move stat producers into PREPROCESSING

`subworkflows/local/preprocessing.nf` — append at the bottom, before `emit:`:

```groovy
// ── Statistics & MultiQC adapters ─────────────────────────────────────────
// AGAT does not handle .gff3.gz transparently — gunzip raw source GFFs first.
GUNZIP_RAW_SOURCE_GFF(
    ch_sources.map { meta, _fa, gff3 -> tuple(meta + [kind: 'raw'], gff3) }
)

ch_agat_in = GUNZIP_RAW_SOURCE_GFF.out.gunzip
    .mix(FILTER_TRANSCRIPT.out.gff3 .map { m, g -> tuple(m + [kind: 'filtered'], g) })
    .mix(RELOCATE_LOCI    .out.gff3 .map { m, g -> tuple(m + [kind: 'decoy'],    g) })

AGAT_SPSTATISTICS(ch_agat_in)
AGAT_TO_MQC(AGAT_SPSTATISTICS.out.stats_yaml)

SEQKIT_STATS(RENAME_FASTA_HEADERS.out.fasta)
SEQKIT_TO_MQC(SEQKIT_STATS.out.stats)

// Emit tuple(meta, path) so downstream consumers can trace provenance.
ch_mqc_files = AGAT_TO_MQC.out.summary
    .mix(AGAT_TO_MQC.out.full)
    .mix(SEQKIT_TO_MQC.out.mqc)
```

Add `mqc_files = ch_mqc_files` to `emit:`. New includes at the top of the
file:

```groovy
include { GUNZIP as GUNZIP_RAW_SOURCE_GFF } from '../../modules/nf-core/gunzip/main'
include { AGAT_SPSTATISTICS               } from '../../modules/local/agat_spstatistics'
include { AGAT_TO_MQC                     } from '../../modules/local/agat_to_mqc'
include { SEQKIT_STATS                    } from '../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_TO_MQC                   } from '../../modules/local/seqkit_to_mqc'
```

Drop the same includes from the soon-to-be-deleted `stats.nf`.

## 4. Move target-side stats into BENCHMARKING

`subworkflows/local/benchmarking.nf` — add AGAT on the target raw + filtered
GFFs, mix gffcompare stats:

```groovy
include { AGAT_SPSTATISTICS as AGAT_TARGET } from '../../modules/local/agat_spstatistics'
include { AGAT_TO_MQC       as AGAT_TARGET_TO_MQC } from '../../modules/local/agat_to_mqc'

// (existing GUNZIP_TARGET_GFF / FILTER_TARGET / GFFCOMPARE wiring)

// Target raw + filtered AGAT — characterises the gffcompare reference.
ch_target_agat_in = GUNZIP_TARGET_GFF.out.gunzip
    .map  { m, g -> tuple(m + [kind: 'raw'], g) }
    .mix(FILTER_TARGET.out.gff3.map { m, g -> tuple(m + [kind: 'filtered'], g) })

AGAT_TARGET(ch_target_agat_in)
AGAT_TARGET_TO_MQC(AGAT_TARGET.out.stats_yaml)

ch_mqc_files = GFFCOMPARE.out.stats
    .mix(AGAT_TARGET_TO_MQC.out.summary)
    .mix(AGAT_TARGET_TO_MQC.out.full)

emit:
    stats     = GFFCOMPARE.out.stats   // kept for any direct consumer
    mqc_files = ch_mqc_files
```

This is the new home for target raw GFF stats — they belong with the
subworkflow that consumes the target GFF.

Because both PREPROCESSING and BENCHMARKING emit AGAT files using the same
`agat_summary` / `agat_full` section IDs, the MultiQC report shows **one
combined table** with one row per (sample, kind) — source and target side
by side, distinguished by sample name (e.g. `Lycaena_phlaeas_282391.raw`
vs. `Lycaena_alciphron_282377.filtered`).

## 5. Split `assets/multiqc_config.yml` into two files

Create `assets/multiqc/` and split the monolithic config into two files
that map onto MULTIQC's two config slots:

- **`main.yml`** — top-level report metadata. Goes into MULTIQC's
  `multiqc_config` slot:
  - `report_title`, `report_comment`
  - `extra_fn_clean_exts` (union of all stages — kept here so we don't
    chase the list across multiple files)
  - `custom_content.order`
  - `top_modules`

- **`sections.yml`** — section/data definitions. Goes into MULTIQC's
  `extra_multiqc_config` slot:
  - `custom_data.agat_summary`, `custom_data.agat_full`,
    `custom_data.seqkit_stats`
  - `sp.agat_summary`, `sp.agat_full`, `sp.seqkit_stats`
  - Any future section IDs for new tools (samtools custom blocks, etc.)
    land here too.

This split survives growth: as new tools come online they add a
`custom_data` block + `sp` pattern to `sections.yml` without touching
`main.yml` (low-risk edits). Top-level layout knobs stay in `main.yml`
(rare, high-impact edits).

`nextflow.config` — retire the single `multiqc_config` param, introduce
two:

```groovy
params {
    // ...
    multiqc_main_config     = "${projectDir}/assets/multiqc/main.yml"
    multiqc_sections_config = "${projectDir}/assets/multiqc/sections.yml"
}
```

(Remove the old `multiqc_config` line.)

Note on AGAT table sharing: both PREPROCESSING and BENCHMARKING emit
`*_summary_mqc.tsv` + `*_full_mqc.tsv` files. Because they use the same
section IDs (`agat_summary`, `agat_full`), MultiQC merges them into one
combined table — source and target rows side by side. Only one
`custom_data.agat_*` block is needed in `sections.yml`.

## 6. Top-level wiring

`workflows/fomo.nf`:

```groovy
include { PREPROCESSING } from '../subworkflows/local/preprocessing'
include { PROJECTION    } from '../subworkflows/local/projection'
include { BENCHMARKING  } from '../subworkflows/local/benchmarking'
include { REPORTING     } from '../subworkflows/local/reporting'

workflow FOMO {
    // ...samplesheet branch as today...

    PREPROCESSING(ch_input.source)

    PROJECTION(
        ch_input.target,
        PREPROCESSING.out.spliced_fasta,
        PREPROCESSING.out.decoy_spliced_fasta
    )

    BENCHMARKING(ch_input.target, PROJECTION.out.gff3)

    REPORTING(
        PREPROCESSING.out.mqc_files
            .mix(BENCHMARKING.out.mqc_files)
            // .mix(PROJECTION.out.mqc_files) — add when PROJECTION grows stats
    )
}
```

STATS's 7-channel signature collapses to REPORTING's single `mqc_files`
input. Adding a new stat producer becomes one `.mix(...)` line — no
REPORTING signature change.

## 7. Module-config selector renames

`conf/modules.config`:

- `FOMO:STATS:GUNZIP_RAW_GFF`     → `FOMO:PREPROCESSING:GUNZIP_RAW_SOURCE_GFF` *and* (new) `FOMO:BENCHMARKING:GUNZIP_TARGET_GFF` already exists, no rename there.
- `FOMO:STATS:AGAT_SPSTATISTICS`  → `FOMO:PREPROCESSING:AGAT_SPSTATISTICS` + `FOMO:BENCHMARKING:AGAT_TARGET`.
- `FOMO:STATS:AGAT_TO_MQC`        → `FOMO:PREPROCESSING:AGAT_TO_MQC` + `FOMO:BENCHMARKING:AGAT_TARGET_TO_MQC`.
- `FOMO:STATS:SEQKIT_STATS`       → `FOMO:PREPROCESSING:SEQKIT_STATS`.
- `FOMO:STATS:SEQKIT_TO_MQC`      → `FOMO:PREPROCESSING:SEQKIT_TO_MQC`.
- New: `FOMO:REPORTING:MULTIQC`.

`ext.prefix` closures unchanged — they're driven by `meta`, not by selector
path.

## 8. CLAUDE.md update

Add a short "Reporting conventions" section:

```
## Reporting conventions

Every subworkflow that emits statistics:
- Runs its stat producers (AGAT, SeqKit, samtools stats, gffcompare, ...)
  and any required adapter modules (modules named *_TO_MQC that emit
  `*_mqc.tsv` files with embedded `# pconfig:` headers).
- Emits a `mqc_files` channel of shape `tuple(meta, path)` — same shape as
  the rest of the pipeline so the meta is available for tracing / logging.
- Subworkflows with no stats emit `Channel.empty()` as `mqc_files`.

The top-level workflow mixes `mqc_files` across subworkflows and passes
the union to REPORTING. REPORTING strips meta, collects the paths, and is
the only place that calls MULTIQC.

MultiQC config lives in two files under `assets/multiqc/`:
- `main.yml` — top-level layout (title, comment, `extra_fn_clean_exts`,
  ordering). Fed via MULTIQC's `multiqc_config` slot.
- `sections.yml` — `custom_data` blocks + `sp` patterns for every custom
  section. Fed via MULTIQC's `extra_multiqc_config` slot.

When a new tool produces statistics, the new MultiQC section's
`custom_data` + `sp` entries go into `sections.yml`. New filename suffixes
to clean are added to `extra_fn_clean_exts` in `main.yml`.
```

## Verification

1. `nextflow run . -profile test,docker -resume`. Expect:
   - All non-STATS tasks cached.
   - PREPROCESSING re-runs only its new AGAT/SeqKit additions (which are
     the same processes that ran before under STATS — Nextflow may cache
     them across the rename if work dir hashes match).
   - BENCHMARKING re-runs the new target AGAT tasks (3 total: 1 raw + 2
     filtered).
   - REPORTING runs MULTIQC_MERGE_CONFIG + MULTIQC once.
2. Compare `results/multiqc/multiqc_report.html` against the pre-refactor
   report (kept in git history or stashed):
   - Same number of `agat_summary` rows (source raw + filtered + decoy +
     target raw + target filtered).
   - Same number of `agat_full` rows.
   - Same `seqkit_stats` rows.
   - Same gffcompare section.
3. `cat results/multiqc/multiqc_data/multiqc_software_versions.txt` —
   confirm same tool set as before.
4. Diff `assets/multiqc/main.yml` + `sections.yml` against the
   pre-refactor `multiqc_config.yml` — set difference should be only
   whitespace + comments.
5. **Negative test**: temporarily drop `BENCHMARKING.out.mqc_files` from
   the `.mix()`. Confirm REPORTING still runs (with fewer files), and the
   resulting report lacks the gffcompare section but is otherwise intact.
   Confirms the new convention is independent per producer.

## Phasing (optional)

If this lands as one PR, the diff is large but mechanical. If staged:

- **Phase A:** add `mqc_files` emit to PREPROCESSING (still keep STATS for
  now — STATS just consumes the new channel + the existing target+raw+
  benchmarking inputs). No reporting changes user-visible. Safe.
- **Phase B:** add `mqc_files` emit to BENCHMARKING; collapse STATS's
  benchmarking-side knowledge.
- **Phase C:** rename STATS → REPORTING; drop the now-empty `take:`
  channels.
- **Phase D:** split `multiqc_config.yml` into `main.yml` + `sections.yml`;
  switch REPORTING to use both MULTIQC config slots.

I'd land A+B+C as one PR (refactor), D as a follow-up PR (config split,
no functional change).

## Out of scope (defer until later)

- **Per-stage MultiQC reports.** The two-config approach still yields one
  report. Splitting into `preprocessing.html` / `projection.html` /
  `benchmarking.html` is a follow-up — REPORTING would call MULTIQC
  multiple times with stage-scoped file subsets, and the YAML split would
  go from main+sections to one config per stage.
- **Per-stage YAML split.** Today we have 2 files (main + sections);
  splitting `sections.yml` further into per-stage files is the natural
  next step once a single sections file gets unwieldy — paired with the
  per-stage report change above.
- **Shared adapter library (`bin/mqc_lib.py`).** Two adapters today don't
  justify a shared library. Revisit at 4+.
- **`extra_fn_clean_exts` shrinkage.** Stays as a union list in `main.yml`
  until per-stage reports naturally scope it down.
- **Cross-stage headline summary report.** Useful but separate.
- **`stats.nf:16` stale comment.** Gets deleted as part of the rename, no
  separate fix needed.
