# 15 — Single mapping: one all-sources alignment per target

## Problem

FOMO aligns **every source separately onto every target**: `F·S·D·T` minimap2 runs, each
re-loading the target `.mmi`, each followed by its own `BAM_TO_GFF`, its own projected-lncRNA
TD2 pass and its own `SAMTOOLS_STATS`. The per-target consensus is then assembled by handing
`gffcompare` the *N* per-source GFF3s (`GFFCOMPARE_COMBINE` in `subworkflows/local/projection.nf`),
and `top3` the same way from a filtered subset (`subworkflows/local/consensus_top.nf`).

Two costs follow. The index is re-read S times per target, and there is **no raw aggregate
annotation** — the only all-source model in the output is the gffcompare-collapsed one, so
"how good is the union of everything before collapsing?" cannot be answered.

## What changes

Invert the flow. Merge the source spliced FASTAs **once** per feature class, run **one**
minimap2 per target, convert that single BAM to a single `allModels.raw.gff3`, and *derive*
everything else from it by splitting on the source species already encoded in the read name.

The species is already carried end-to-end and needs no new plumbing:
`modules/local/rename_fasta_headers.nf` writes `>tid|type|species`, minimap2 puts that in the
BAM QNAME, and `bin/bam_to_gff.sh` stamps it onto every `gene`/`transcript` row as
`source=<species>`. Splitting is a text operation on an attribute that already exists.

```
PREPROCESSING (unchanged, per source)
   └─ *.spliced.renamed.fasta          headers: >tid|type|species

PROJECTION
   MERGE_SOURCE_FASTA        F·D tasks, target-independent
      └─ allsources.<gtype>.fasta
   MINIMAP2_INDEX            T tasks (unchanged)
   MINIMAP2_ALIGN            F·D·T tasks   ← was F·S·D·T
      └─ <target>.allModels.<ft>[.decoy].raw.bam(.bai)             ★ published
           ├─ SAMTOOLS_STATS / SAMTOOLS_TO_MQC
           └─ BAM_TO_GFF  (intermediate, unpublished)
                └─ EXTRACT_PROJECTED_LNC → TD2_NONCODING_PROJ → FILTER_ALLMODELS
                     └─ <target>.allModels.<ft>[.decoy].raw.gff3   ★ published
                          ├─ SPLIT_GFF_BY_SOURCE  (new script)
                          │     └─ <target>.from_<src>.<ft>[.decoy].projected.gff3  ★ published
                          └─ GFFCOMPARE_COMBINE → COMBINED_GTF_TO_GFF
                                └─ <target>.allModels.<ft>[.decoy].collapsed.gff3   ★ published

BENCHMARKING   gffcompare vs target reference for: each source, allModels.raw, allModels.collapsed

CONSENSUS_TOP
   SELECT_TOP_SOURCES  (per-source lncRNA F1, self excluded — unchanged)
      └─ SUBSET_GFF_BY_SOURCE on allModels.raw
           └─ <target>.top3.<ft>.raw.gff3                          ★ published
                └─ GFFCOMPARE_COMBINE_TOP → COMBINED_TOP_GTF_TO_GFF
                     └─ <target>.top3.<ft>.collapsed.gff3          ★ published
   GFFCOMPARE_TOP on top3.raw AND top3.collapsed
```

Accuracy values per annotated target: **S** per-source + `allModels.raw` +
`allModels.collapsed` + `top3.raw` + `top3.collapsed`.

`gffcompare` in combine mode with a **single** input file is already exercised today — the
test run's `Lycaena_alciphron` consensus has exactly one qualifying source and still yields
`.combined.gtf` — so collapsing one merged GFF3 needs no new mechanism.

### Decisions

| Question | Decision | Rationale |
|---|---|---|
| Split level | **GFF**, on the `source=` attribute | No per-source BAMs to write; `bam_to_gff.sh` is record-local so the split is a pure line filter |
| Self-projection in `allModels` | **Included** | `allModels` means literally all sources |
| Self-projection in top3 ranking | **Excluded** (unchanged from plan 14) | Self always wins on F1; top3 would restate the target's own annotation |
| Projected-lncRNA TD2 filter | **Once on `allModels`**, split runs on the *filtered* GFF | S·T TD2 runs → D·T |
| `--include_mrna` / `--include_decoy` | **One merged FASTA per feature class** (`lnc_RNA`, `mRNA`, `lnc_RNA_decoy`, `mRNA_decoy`) | Keeps the gene-type separation every downstream grouping already assumes; no class demux in the split script |

### Cardinality

S sources, T targets, T_g ≤ T with a `gff3`, F feature types, D = 2 with `--include_decoy`:

| Stage | Before | After |
|---|---|---|
| PREPROCESSING | `F·S` | `F·S` (unchanged) |
| MERGE_SOURCE_FASTA | — | `F·D` |
| MINIMAP2_INDEX | `T` | `T` |
| MINIMAP2_ALIGN / SAMTOOLS_STATS / BAM_TO_GFF | `F·S·D·T` | **`F·D·T`** |
| projected TD2 (`EXTRACT` + `TD2_PREDICT` + `TD2_CODING_FILTER`) | `S·T` | **`T`** (lncRNA only, and now one pooled run per target) |
| FILTER_ALLMODELS (was FILTER_PROJ_LNC_GFF) | `S·T` | `F·D·T` (every class, non-lncRNA with an empty drop list) |
| SPLIT_GFF_BY_SOURCE | — | `F·D·T` (emitting `F·S·D·T` files) |
| GFFCOMPARE_COMBINE | `F·T·D` (N inputs) | `F·T·D` (1 input) |
| BENCHMARKING GFFCOMPARE | `(F·S·D + F·D)·T_g` | `(F·S·D + 2·F·D)·T_g` |
| GFFCOMPARE_TOP | `F·T_g` | `2·F·T_g` |

### Two accepted consequences

1. **Per-source published GFF3s stop being byte-identical to `results/`.** The prompt asked
   for byte-identity; the chosen TD2 placement (once on `allModels`, single split) means the
   split runs *after* the coding filter, so per-source files lose exactly the lncRNA
   transcripts TD2 calls coding. `mRNA` and decoy tracks see no TD2 and **do** stay
   byte-identical. Pooling also shifts a handful of coding calls, because TD2.Predict trains
   on its input set.
2. **MultiQC loses per-source alignment rows.** `SAMTOOLS_STATS` now runs on the merged BAM:
   one row per target per feature class. Per-source *model* counts are unaffected —
   `GFF_STATS_PROJECTED` still runs on every per-source GFF3.

### Risk to watch: duplicate transcript IDs

Per-source GFF3s were separate files, so a transcript id shared by two species never
collided. In `allModels.raw.gff3` they land in one file, where a duplicate `ID=` would make
`gffread -w` and the `gffcompare` collapse silently misbehave. Ensembl ids are species-scoped
so this is unlikely, but `bin/gff_by_source.py` **asserts uniqueness and exits 1** rather than
letting it pass.

## Implementation phases

Split into three, each independently runnable and verifiable. Phase A is deliberately
**behaviour-preserving downstream**: it reconstitutes the per-source GFF channel with the same
shape and meta, so benchmarking and consensus keep working untouched and any regression is
isolated to the alignment change.

| Phase | Plan | Delivers | Status |
|---|---|---|---|
| A | [`15a_single_mapping_align.md`](15a_single_mapping_align.md) | Merged FASTA, one alignment per target, `allModels.raw.gff3`, split back to per-source. Downstream unchanged. | ✅ done, checkpoint passed 2026-08-11 |
| B | [`15b_single_mapping_aggregates.md`](15b_single_mapping_aggregates.md) | Collapse from `allModels`, `top3.raw` by subset, and the four aggregate accuracy values. | ✅ done, checkpoint passed 2026-08-12 |
| C | [`15c_single_mapping_outputs.md`](15c_single_mapping_outputs.md) | Output layout, reporting, resources, docs. | ✅ done, checkpoint passed 2026-08-12 |

Each phase ends with a **Checkpoint** that re-checks the code against both that phase's scope
and this document. Do not start the next phase until its checkpoint passes.

## Baseline

`results/` is the output of `nextflow run main.nf -profile test,crg` on the current code.
**Copy it to `results_baseline/` before the first new run** — every phase diffs against it.
