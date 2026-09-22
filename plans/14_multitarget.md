# 14 — Multi-target: `role: both` and all-vs-all in one run

## Problem

FOMO annotated **one** target from **many** sources. An all-vs-all experiment therefore
meant N separate `nextflow run` invocations with N hand-generated samplesheets — 20
species, 20 runs, 20 `--outdir`s. That repeated all source-side preprocessing N times
(it is entirely target-independent) and scattered results across N output trees.

## What changed

The samplesheet now describes the whole experiment:

- `role` gains **`both`** — donor *and* target.
- `gff3` is **optional for a pure `target`** (an un-annotated assembly can be annotated
  but not benchmarked); still required for `source`/`both`, enforced by an
  `allOf`/`if`-`then` in `assets/schema_input.json`.

One run now produces every target × source projection, with preprocessing done once per
species.

### Decisions

| Question | Decision | Where |
|---|---|---|
| Self-pairs (X→X) | **Run** them — an Sn/Pr ceiling control — but **exclude from the `combined` consensus and the top-N ranking pool**, else self always ranks #1 on F1 and the consensus degenerates into a copy of the target's own annotation | `projection.nf` (combine input filter), `consensus_top.nf` (ranking + pool filters) |
| Target with no `gff3` | Per-source projections **and** the all-source `combined` consensus; no gffcompare, no target GFF stats, no top-N | `benchmarking.nf` filters `ch_target` on a truthy `gff3`; the inner join drops the rest |
| MultiQC | **One report per target**; target-agnostic source-side sections broadcast into each | `reporting.nf` |
| Output layout | `results/targets/<target>/<process>/…` | `nextflow.config` publishDir |

### Cardinality

S sources, T targets, T_g ≤ T with a `gff3`, F feature types, D = 2 with `--include_decoy`:

| Stage | Tasks |
|---|---|
| PREPROCESSING | `F·S` filter, `S` decoy chains — **scales with S, not S·T** |
| MINIMAP2_INDEX | `T` |
| MINIMAP2_ALIGN / BAM_TO_GFF | `F·S·T·D` |
| GFFCOMPARE_COMBINE | `F·T·D` |
| BENCHMARKING GFFCOMPARE | `(F·S·D + F·D)·T_g` |
| CONSENSUS_TOP | one `SELECT_TOP_SOURCES` + `F` combines per `T_g` |
| MULTIQC | `T` |

## Implementation notes (the non-obvious parts)

**Do not nest the subworkflows.** The obvious implementation — a `PER_TARGET` subworkflow
looped over targets — breaks ~30 `ext.prefix` / `ext.sample_name` closures **silently**:
every selector in `conf/modules.config` is a fully-qualified process path and `withName`
matching is a regex *find*, so an extra level makes the substring stop matching and modules
fall back to `${meta.id}` (the **source** id). The target dimension lives in the channels.

**Two `.filter`s, not `.branch`** (`workflows/fomo.nf`). `branch` routes each item to
exactly one output, so a `both` row could never reach both channels. The samplesheet is
also read into a plain List first, so three invariants the JSON schema cannot express
(unique species, ≥1 target, ≥1 source) abort at launch with a clear message.

**Every projection ⋈ reference join keys on `[target_id, feature_type]`**
(`benchmarking.nf`, `consensus_top.nf`). Keyed on `feature_type` alone, each projection is
also compared against every *other* target's annotation — and since `ext.prefix` is built
from the **query** meta, all T of those tasks emit the same filename into the same
published directory. Wrong numbers under a plausible-looking name.

**`SELECT_TOP_SOURCES` is grouped per target.** `bin/select_top_sources.py` keys its score
map by source alone, so pooling several targets' stats into one task would have them
overwrite each other. Grouping in the channel leaves the script unchanged; self-exclusion
is done in the channel too, for the same reason.

**`MINIMAP2_ALIGN`'s reference is now a per-task `multiMap` branch**, not a broadcast value
channel. The two `.first()` calls in `projection.nf` are gone — the second one (target
FASTA for `EXTRACT_PROJECTED_LNC`) was the dangerous one: with T > 1 it would extract
projected transcript sequences out of the wrong genome and hand TD2 nonsense to score.

**REPORTING routes by the *presence* of `meta.target_id`** — set → that target's report
only; absent → broadcast to every report. Shared files are fanned out **one file per item**
(`flatMap`, not `map`) against the target ids. Two reasons, and skipping either fails at
runtime: `GFF_STATS_TO_MQC` emits a *list* of two TSVs per item, and `combine` SPREADS a
List-valued item across the output tuple — so `[meta, [a, b]].combine(ids)` yields
`[a, b, tid]` and the next closure is invoked with three arguments
(`Invalid method invocation 'call' with arguments: [...] on _closure8 type`). Same reason
you cannot `.collect()` the shared files and combine that list in.

**Target-side GFF stats gained a `.target` / `_target` segment** in `conf/modules.config`.
With `role: both` the same species writes a source row and a target row into the same
`gffstats_input` table, and they are genuinely different data — the source-side lncRNA is
TD2 coding-filtered, the target reference deliberately is not. Without the segment both
sides collide on filename and row key, and one silently clobbers the other.

**`TD2_PREDICT` needed a second guard** (`modules/local/td2_predict.nf`), surfaced by
enabling self-projection. It already skipped an empty input FASTA; it now also skips when
`TD2.LongOrfs` extracts **zero complete ORFs** from a non-empty one. That is the *norm* for
a self-projection — its lncRNA are precisely the transcripts the upstream PREPROCESSING TD2
pass already cleared of complete ORFs — and it left PSAURON with nothing to score, no
`psauron_score.csv`, and `TD2.Predict` dying in `pandas.read_csv`. Pre-existing latent bug
(any input with no complete ORF hits it); self-pairs made it near-certain.

## Test profiles

- `-profile test` → `assets/samplesheet.csv`: 2 rows `both` (→ self-pairs), 1 pure `target`
  with an **empty** `gff3`, 2 pure `source`. T = 3, S = 4 → 12 pairs.
- `-profile test_single` → `assets/samplesheet_single.csv`: the pre-multitarget layout
  (1 target, 4 sources, gff3 everywhere), so the T = 1 path stays covered.

Both run **lncRNA only** (F = 1, D = 1) — the mRNA positive control and the decoy FDR
baseline are development instruments, and enabling both quadrupled the task count without
testing anything the lncRNA path does not already cover. Pass `--include_mrna` /
`--include_decoy` explicitly to exercise those branches.

On the cluster use `-profile test,crg` — in that order, since both profiles set
`process.resourceLimits` and the later one wins.
