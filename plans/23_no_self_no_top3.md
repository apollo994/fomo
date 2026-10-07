# 23 — Remove self-projection, the top-N consensus and allModels_collapsed

## Status: implemented (2026-10-07) on branch `no-self-no-top3`, off `annotated-targets` @ `99f648b` — see **Results**

## Decisions (2026-10-07)

| Question | Decision |
|---|---|
| Self-projection (X→X for a `both` species) | **Removed entirely**, for every track (lncRNA, mRNA, decoys). A species' own sequences are never aligned to its own genome. |
| Top-N consensus (`top3`, source ranking) | **Removed entirely.** The ranking is biased by whether a reference exists and by its quality. The cross-species-supported curated set (plans/21) replaces it as the "consensus". |
| `allModels_collapsed` (gffcompare-combine of allModels) | **Removed** (process, published file, stats, benchmark). |
| What is benchmarked (gffcompare) | **Only `allModels` raw and each per-source model set** (`from_<source>`), per gene type, including decoys as today. Curated / supported models are **not** benchmarked. |
| Benchmark reference | **Unchanged**: the target's filtered ★ annotation per gene type (spliced, longest isoform, lncRNA TD2-noncoding; plans/20), which exists for `both` rows only. |
| Which targets are benchmarked | **Only `both` rows** (plans/22 unchanged: a pure target's gff3 is the curation reference only). |
| Run-summary "added candidates" (was allModels_collapsed / top3_collapsed) | **Curated genes** (plans/21) per target and track, which also folds in plans/21's pending phase C. With decoys on, add a curated-decoy rate. *To confirm in review; see Open questions.* |

## Why

- **Self-projection** contributes nothing to the curated deliverable: it never counts as support,
  and its models overlap the reference they came from. Its only job was the Sn/Pr ceiling control.
  It costs a lot:
  - **share of `allModels.lnc_RNA.raw`:** median 57 % on acari, 40 % apiaceae, 20 % solanaceae,
    12 % drosophilidae (measured 2026-10-07 on the existing runs);
  - **inflated `allModels` accuracy:** CLAUDE.md already flags it as "pulled up towards the self
    ceiling";
  - **wasted work:** alignment, TD2, split, stats and benchmarking all process it.
- **Top-N** ranks donors by F1 against the target's own annotation. Targets without a reference, or
  with a poor one, get no or misleading picks, and the "best" donors are those closest to the
  reference's own biases.
- **allModels_collapsed** is a redundancy-collapsed view nobody consumes once top-N and its
  "added" counts are gone. The curated set is the non-redundant deliverable.

## Design

### 1. Self-exclusion at alignment (`modules/nf-core/minimap2/align/main.nf`, FOMO-local branch)

**Where.** `MERGE_SOURCE_FASTA` stays target-independent (F·D tasks), and no new process is
added. The FOMO-local `bam_format && !bam_input` branch of `MINIMAP2_ALIGN` (already a documented
local deviation) gets an optional exclusion step:

```
# meta.exclude_source set → drop query records whose header species (last |-field) == it
awk -v s="${meta.exclude_source}" '/^>/{n=split(substr($1,2),f,"|"); keep=(f[n]!=s)} keep' <reads> > query.fa
minimap2 … ${target} query.fa … -o ${prefix}.unsorted.sam
```

Details:
- **Input.** `zcat -f` first if the merged FASTA is gzipped.
- **Temp file, not a pipe.** The script writes a temporary file and does not pipe into minimap2,
  for the reason the branch already documents: a pipe masks minimap2's own exit code (137 on OOM),
  which the retry strategy keys on. The temp file is removed after alignment.
- **Header grammar.** `>tid|type|species` (RENAME_FASTA_HEADERS); the species is the last
  `|`-field. That field is the chain-of-custody key that `bam_to_gff.sh` turns into `source=`.
- **Who gets the flag.** `projection.nf` sets `exclude_source: tid` on the align meta only when
  the target is also a source (`role == 'both'`). This is decided from `ch_target`'s meta, so pure
  targets skip the copy.
- **Empty query.** When the target's only source is itself, the filtered query is empty. minimap2
  on an empty FASTA must yield a header-only BAM (verify). That is the existing "nothing projected"
  path (CLAUDE.md merge/split contract: `optional: true` split, `join(remainder: true)`), and the
  TD2 zero-ORF guards stay.
- **Container.** The `minimap2_samtools` Wave image must provide `awk` and `zcat`; check with
  `singularity exec`. If one is missing, rebuild the image with coreutils/gawk added (same spec
  otherwise).

Alternatives rejected:
- **Per-target merged FASTAs** (F·D·T merges): this duplicates the largest file T times on disk.
- **Drop self after alignment:** the alignment compute is still spent.

### 2. Remove the top-N consensus
- **Delete:**
  - `subworkflows/local/consensus_top.nf`;
  - `modules/local/select_top_sources.nf`, `bin/select_top_sources.py`;
  - the `--keep` subset mode of `bin/gff_by_source.py` and `modules/local/gff_by_source.nf`.
    It was only used by `SUBSET_GFF_BY_SOURCE`. The split mode stays.
- **`workflows/fomo.nf`:**
  - remove the `CONSENSUS_TOP` call;
  - `REPORTING` takes `BENCHMARKING.out.stats` only;
  - `RUN_SUMMARY` loses its `top_sources` input.
- **Params:** remove `top_consensus_n` from `nextflow.config`, `nextflow_schema.json` and the
  launch summary.
- **`conf/modules.config`:**
  - remove the `CONSENSUS_TOP` block and the `select_top_sources/` publish rule;
  - drop `_TOP` from the alternations `(GFFCOMPARE|GFFCOMPARE_BATCH|GFFCOMPARE_TOP)$` and
    `(…|SUBSET_GFF_BY_SOURCE|COMBINED_TOP_GTF_TO_GFF)$`;
  - update the docker-profile `GFF_STATS` selector in `nextflow.config`
    (`'.*:GFF_STATS(_TARGET|_PROJECTED(_BATCH)?|_TOP)?$'` → without `_TOP`).
- **`conf/base.config`:** drop any `withName` entries for removed processes.
- **MultiQC:** remove top3-related sections and section-order entries in `main.yml` /
  `sections.yml` and in `run_summary_main.yml` / `run_summary_sections.yml`.

### 3. Remove `allModels_collapsed`
- **`projection.nf`:**
  - remove `GFFCOMPARE_COMBINE` and `COMBINED_GTF_TO_GFF`, and the `allmodels_collapsed` emit;
  - the per-target batch list becomes `[every source's split GFF3] + [allModels raw]`, so drop
    the `ch_collapsed_one` join (`remainder: true` stays for the split side).
- **`modules/local/gff_stats_projected_batch.nf`:** drop the `*.collapsed.gff3` symlink loops.
- **`conf/modules.config`:** remove the `GFFCOMPARE_COMBINE` / `COMBINED_GTF_TO_GFF` blocks and
  their publish alternation entries.
- **`modules/local/gffcompare_batch.nf`:** comments and naming only. It loops over whatever list
  it gets.

### 4. Aggregate-id grammar (the "three places that must agree", CLAUDE.md)
- **`bin/fomo_stats.py`:** `AGGREGATE_IDS = ("allModels_raw",)`; update its docstring, which
  points at `consensus_top.nf`.
- **`benchmarking.nf`:** the re-key closure is unchanged (by length), but its comment lists the
  pseudo-sources.
- **The consensus_top channel filter** goes away with the file, so two places remain.
- **`bin/gffcompare_accuracy_mqc.py`:** only `allModels` raw as the aggregate point; update the
  docstring.

### 5. Run summary (`bin/run_summary_tables.py`, `subworkflows/local/run_summary.nf`)
- **Inputs:** remove `--top-sources` / `_read_top_sources`. Add the curation report files
  (`*.curation.json` from `CURATION.out.report`, or the three `_curation_*_mqc.tsv` already in
  `ch_all_mqc`; prefer the mqc TSVs, which need no new channel).
- **`section_overview`:**
  - drop "… of which self-pairs" / `n_self_pairs` and the top-N param row;
  - "Source × target pairs" becomes S·T − |both|;
  - keep the plans/22 rows.
- **`section_targets`:** columns become
  `before | curated_genes | after | pct_increase | reads_mapped_percent`, per transcript type, with
  `curated_genes` from the curation novel table. The bar plot has a single dataset.
- **`section_species` / `section_donors`:** remove the "times picked in top-N" counts and the
  `run_top_sources_mqc.tsv` table. Keep per-donor accuracy (from the per-source gffcompare rows).
- **`section_aggregates`:** "allModels raw vs best single source" (no top3 / collapsed rows).
- **`section_heatmap`:** unchanged (the source × target matrix now has no diagonal self cells).
- **New, with decoys on:** per target, curated decoy genes / curated real genes, plus the same
  normalised per input model, since `decoy_cap` puts the two on different scales.

### 6. Curation (`bin/curate_models.py`)
- Self exclusion stays as a guard. If any model with `source=<target>` reaches it, log a warning
  ("self-projection present, input predates plans/23") and exclude it as today. `--keep-self`
  stays for standalone use on old runs.
- Docstring: drop the top3 mention.

### 7. Docs
- **`CLAUDE.md`:**
  - rewrite the "Self-pairs" paragraph: self-projection is removed, and the ceiling control is
    gone with it;
  - the "allModels means literally every source" statement gets "except the target itself";
  - subworkflow map: no `CONSENSUS_TOP`;
  - PROJECTION emits;
  - task-count paragraph: no `GFFCOMPARE_COMBINE`, `GFFCOMPARE_TOP`, `2·F·T_g`; benchmark
    comparisons become `(F·(S−1)·D + F·D)·T_g` for `both` targets;
  - meta-map contract `id` values become `allModels | allModels_raw`;
  - merge/split contract: drop subset mode and the `transpose()` single-file note;
  - the TD2 zero-ORF paragraph: the guard now fires when nothing but the target itself would have
    projected;
  - output tree: no `top3.*`, no `allModels.*.collapsed.gff3`, no `select_top_sources/`;
  - "Reporting conventions" mentions of top3.
- **`README.md`:** pipeline stages (no consensus step; curation is the consensus); outputs.
- **`conf/test.config`:** the comment about self-pairs.
- **Older plans:** plans/09, 15, 18 and 20 are not edited (history). This plan supersedes their
  top-N and self-pair parts.

## Task counts after this plan (per run)

The `F·D·T` alignments stay. Removed:
- `F·D·T` `GFFCOMPARE_COMBINE`;
- `F·D·T` `COMBINED_GTF_TO_GFF`;
- `T_g` `SELECT_TOP_SOURCES`;
- `F·T_g` `SUBSET_GFF_BY_SOURCE`;
- `F·T_g` `GFFCOMPARE_COMBINE_TOP` + `F·T_g` `COMBINED_TOP_GTF_TO_GFF`;
- `2·F·T_g` `GFFCOMPARE_TOP` + top GFF stats.

Each `both` target's alignment, TD2, split and stats lose the self share of their input.

## Verification

1. **`-profile test,crg`** (the plans/22 samplesheet: alciphron `both`, thersamon `source`,
   hippothoe `target` + gff3):
   - no `from_Lycaena_alciphron_282377` file under alciphron;
   - no `source=Lycaena_alciphron_282377` in its `allModels.*.raw.gff3`;
   - no `top3.*`, `*.collapsed.gff3`, `select_top_sources/`;
   - alciphron's gffcompare has per-source (thersamon) + `allModels_raw` rows only;
   - hippothoe is not benchmarked;
   - reports and the run summary render;
   - run-summary "curated genes" equal the curation novel tables;
   - no `targets/null/`.
2. **Same with `--include_mrna --include_decoy`:** decoys excluded from self too (no
   `source=alciphron` in its decoy allModels); curated-decoy rate present.
3. **Empty-query edge:** a sheet where a `both` species is the only donor (alciphron `both` +
   hippothoe `target`). Alciphron's query is empty after exclusion, so it gets an honest
   0-model result and the run completes; hippothoe still gets alciphron's models.
4. **Regression against a pre-23 run** of the same sheet:
   - per-source gffcompare stats for non-self pairs are identical or explained. minimap2 maps
     each query independently, but the projected-lncRNA TD2 pass is pooled per target, and
     removing self changes its training set. Expect at most a few coding calls to differ
     (CLAUDE.md, plans/15);
   - curated outputs are identical up to that same TD2 effect.
5. `nextflow lint .`: no new errors (the 7 known on `main`); the `test_single` include error is
   unrelated.
6. **Container check:** `awk` and `zcat` are present in the minimap2_samtools image.

## Open questions (defaults in **bold**)

1. **Run-summary "added" metric** = **curated genes per target and track** (plans/21 phase C
   folded in). *Implemented as the default.*
2. **Per-donor reporting without top-N:** **keep the per-donor accuracy (F1 per source × target)
   heatmap and the donor table without "picked" counts.** *Implemented as the default.* Alternatively, add "curated genes
   supported by this donor" (how often a donor's models appear in `supporting_models`), which
   would be a reference-free measure of donor usefulness.

## Results (2026-10-07)

**Implementation notes vs the design:**
- **Self-exclusion** is in the FOMO-local branch of `modules/nf-core/minimap2/align`:
  - a mawk filter on the last `|`-field of each header writes `<prefix>.query.fa`, which is
    removed after alignment;
  - it is triggered by `meta.exclude_self` (set in `projection.nf` from the target's
    `role == 'both'`);
  - the minimap2_samtools image has mawk 1.3.4 and zcat, so no image rebuild was needed;
  - minimap2 on an empty query exits 0 and writes a header-only BAM (verified in the image).
- **"Added" in the run summary is counted in GENES:** `genes_before` is the reference's
  `n_genes` for that type; `curated_genes` is one representative per gene. With decoys,
  `est_curated_fdr_pct` = (decoy curated genes / decoy multi-exon input models) /
  (real curated genes / real multi-exon input models) × 100.
- **Donor table:** `run_donor_ranking_mqc.tsv` keeps its name but holds accuracy only; the
  picks bargraph and the top-sources table are gone.
- **`curate_models.py`** warns (instead of the old "note") only if self models still reach it.

**`-profile test,crg`, default sheet** (alciphron `both`, hippothoe `target` + gff3, thersamon
`source`), not resumed: 65 tasks succeeded.
- **No self-projection:** alciphron's `allModels.lnc_RNA.raw.gff3` has 0 `source=alciphron`
  models (60 transcripts, was 115 with 55 self). There is no `from_Lycaena_alciphron_282377`
  file.
- **Benchmarking:** alciphron's gffcompare has exactly `from_thersamon` + `from_allModels_raw`.
- **Removed outputs:** no `top3.*`, `*.collapsed.gff3`, `select_top_sources/` or
  `*top_sources*`, and no such sections in any report.
- **Regression vs the plans/22 run of the same sheet:**
  - alciphron ← thersamon projected GFF3 is byte-identical, with Transcript Sn/Pr 28.8 / 28.3
    on both runs;
  - alciphron `allModels_raw` went from 91.5 / 47.0 (self-inflated) to 28.8 / 28.3, i.e.
    exactly its one real donor;
  - hippothoe's curated + merged GFF3s and curation report are byte-identical;
  - alciphron's curation report differs only in `n_models` (115 → 60) and `self_excluded`
    (55 → 0).
- **Run summary:**
  - "Source × target pairs (no self-projection)" = 3;
  - targets table in genes: hippothoe lncRNA 64 + 3 curated, alciphron 62 + 0;
  - every section renders once; no duplicated sections in the per-target reports.

**`--include_mrna --include_decoy`** (resumed): 111 new + 60 cached tasks succeeded.
- 0 self models on all four alciphron tracks (lncRNA 60, lncRNA decoy 40, mRNA 830, mRNA decoy
  755 transcripts).
- gffcompare: thersamon + allModels_raw per track.
- Curated decoy genes 0. Hippothoe mRNA +1 curated gene.

**Empty-query edge case** (`alciphron both` + `hippothoe target`, alciphron is the only donor):
- alciphron's query is empty after exclusion, its `allModels.lnc_RNA.raw.gff3` has 0
  transcripts, and there are no split files;
- benchmarking ran on the empty `allModels_raw` (the honest 0-model result);
- curation reports 0 models / 0 genes;
- hippothoe still receives alciphron's models. The run completed.

**`bin/gff_by_source.py`** (split-only now) reproduces all four per-source files of the plans/22
run byte-for-byte.

**`nextflow lint .`:** the same 7 errors as before, none new (41 clean files, 2 fewer:
`consensus_top.nf` and `select_top_sources.nf` were deleted).
