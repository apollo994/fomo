# 15a — Phase A: merge, align once, split back

Sub-plan of [`15_single_mapping.md`](15_single_mapping.md). Read that first.

## Goal

Collapse `F·S·D·T` minimap2 runs into `F·D·T`, produce `allModels.raw.gff3`, and
**reconstitute the per-source GFF3 channel with an identical shape and meta** so that
`BENCHMARKING`, `CONSENSUS_TOP` and `REPORTING` keep working with no edits. Any regression
after this phase is attributable to the alignment change alone.

Explicitly **out of scope** here: `GFFCOMPARE_COMBINE` still consumes the N per-source GFF3s
(self excluded) exactly as today; `CONSENSUS_TOP` is untouched; publish layout is untouched.

## Changes

### New — `bin/gff_by_source.py`

One script, two modes. Reads a GFF3, keys records on the `source=` attribute that
`bin/bam_to_gff.sh` stamps on `gene` and `transcript` rows; `exon` rows inherit through
`Parent`. Records are emitted in input order and `##gff-version 3` is re-written per output.

```
gff_by_source.py --gff IN --outdir out
    --name-template '<target>.from_{source}.<ft>.projected.gff3'   # split mode
    [--keep sources.csv --name '<target>.top3.<ft>.raw.gff3']      # subset mode (Phase B)
```

- Split mode writes one file per distinct `source=` value.
- Subset mode (wired in Phase B, implement now) reads the one-column CSV that
  `bin/select_top_sources.py` already emits (`source` header) and writes a single file.
- **Fails hard on a duplicate transcript `ID=` across the whole input**, with the offending
  id and the two sources in the message.
- Prints per-source record counts to stderr so the `.command.err` is a usable audit trail.

### New — `modules/local/gff_by_source.nf`

```groovy
process GFF_BY_SOURCE {
    tag "${meta.target_id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_single'
    container 'python:3.11'
    input:  tuple val(meta), path(gff3)
    output: tuple val(meta), path("out/*.gff3"), emit: gff3
}
```

Aliased `SPLIT_GFF_BY_SOURCE` in `projection.nf` (and `SUBSET_GFF_BY_SOURCE` in Phase B).
Both consumers `.transpose()` the emitted list, so one module serves both shapes.

### New — `assets/empty_ids.txt`

Zero-byte file. Lets the `mRNA` and decoy tracks pass through the same `FILTER_GFF_BY_ID`
step as lncRNA with an empty drop list, so **one** process produces every
`allModels.*.raw.gff3`. Without it the lncRNA branch would come out of `FILTER_ALLMODELS`
and the other branches out of `BAM_TO_GFF`, which means two publish points and a name
collision on the lncRNA track.

### New — merged-FASTA module

`nf-core modules install cat/cat`, aliased `MERGE_SOURCE_FASTA`. If that install is
inconvenient, a local `modules/local/merge_fasta.nf` (`cat ${fastas} > ${prefix}.fasta`) is
equivalent — what matters is that the input list is sorted by name for determinism.

### Modified — `modules/local/filter_gff_by_id.nf`

Stage the input as `path(gff3, name: 'input/*')` (the trick `modules/local/maybe_gunzip.nf`
already uses) and emit `${prefix}.gff3` instead of `${prefix}.filtered.gff3`, so `ext.prefix`
can own the full deliverable name `<target>.allModels.<ft>.raw`. Append `.filtered` to the
preprocessing-side `FILTER_LNC_GFF` prefix in `conf/modules.config` to keep that (unpublished)
filename stable.

### Modified — `subworkflows/local/projection.nf`

Replace the `ch_spliced_fasta.mix(ch_decoy_spliced).combine(ch_ref)` cartesian fan-out with a
merge-then-fan-out:

```groovy
ch_merged = ch_spliced_fasta.mix(ch_decoy_spliced)
    .map { meta, fa ->
        def gtype = meta.feature_type + (meta.decoy ? '_decoy' : '')
        tuple([id: gtype, feature_type: meta.feature_type, decoy: meta.decoy, gtype: gtype], fa)
    }
    .groupTuple()
    .map { meta, fas -> tuple(meta, fas.sort { it.name }) }   // deterministic concat order

MERGE_SOURCE_FASTA(ch_merged)                                // F·D tasks, target-independent

ch_align = MERGE_SOURCE_FASTA.out.file_out
    .combine(ch_ref)
    .multiMap { meta, reads, tid, mmi, _fa ->
        reads:     tuple(meta + [id: 'allModels', target_id: tid], reads)
        reference: tuple([id: tid], mmi)
    }
```

Keep `multiMap` rather than a broadcast reference, for the reason the existing comment gives:
with T > 1 the reference must be paired with its own reads task-by-task.

Then:

1. `MINIMAP2_ALIGN` / `SAMTOOLS_STATS` / `SAMTOOLS_TO_MQC` / `BAM_TO_GFF` unchanged in
   wiring — they simply now see one item per (target, class).
2. Keep the `branch` → `EXTRACT_PROJECTED_LNC` → `TD2_NONCODING_PROJ` chain, driven by the
   single `allModels` GFF. **Keep its `combine(ch_target_fa, by: 0)` on `target_id`** — each
   projection must be extracted from its own genome.
3. Route **all** classes through `FILTER_ALLMODELS` (alias of `FILTER_GFF_BY_ID`): lncRNA gets
   TD2's `coding_ids`, the rest get `assets/empty_ids.txt`. Output = `ch_allmodels_raw`.
4. `SPLIT_GFF_BY_SOURCE(ch_allmodels_raw)`, then re-key by stripping the *known* prefix and
   suffix — exact string arithmetic, not a regex guess (species names contain `_` and digits):

```groovy
ch_projected_persource = SPLIT_GFF_BY_SOURCE.out.gff3
    .transpose()
    .map { meta, gff ->
        def pre = "${meta.target_id}.from_"
        def suf = ".${meta.feature_type}${meta.decoy ? '.decoy' : ''}.projected.gff3"
        tuple(meta + [id: gff.name[pre.size()..-(suf.size() + 1)]], gff)
    }
```

`ch_projected_persource` then feeds the existing `GFFCOMPARE_COMBINE` grouping, the existing
`GFF_STATS_PROJECTED`, and the existing `emit: gff3` — untouched in this phase.

### Modified — `conf/modules.config`

Fully-qualified selectors only; `withName` is a regex *find*, so never nest the subworkflows.

| Process | `ext.prefix` |
|---|---|
| `FOMO:PROJECTION:MERGE_SOURCE_FASTA` | `allsources.${meta.gtype}` |
| `FOMO:PROJECTION:MINIMAP2_ALIGN` | `${meta.target_id}.allModels.${meta.feature_type}${decoy}.raw` |
| `FOMO:PROJECTION:BAM_TO_GFF` | `${meta.target_id}.allModels.${meta.feature_type}${decoy}.unfiltered` |
| `FOMO:PROJECTION:FILTER_ALLMODELS` | `${meta.target_id}.allModels.${meta.feature_type}${decoy}.raw` |
| `FOMO:PROJECTION:SPLIT_GFF_BY_SOURCE` | `ext.args` supplying `--name-template` |
| `FOMO:PROJECTION:EXTRACT_PROJECTED_LNC` | `${meta.target_id}.allModels.lnc_RNA.raw` |
| `FOMO:PROJECTION:TD2_NONCODING_PROJ:TD2_CODING_FILTER` | prefix/`sample_name` → `${meta.target_id}.allModels` |
| `FOMO:PROJECTION:SAMTOOLS_STATS`, `SAMTOOLS_TO_MQC` | `${meta.target_id}.allModels.${meta.feature_type}${decoy}.raw` |

`GFF_STATS_PROJECTED*` closures keep `from_${meta.id}` — after the split, `meta.id` is the
source species again, so the projection stats table is unchanged.

Add `FILTER_ALLMODELS` to the publishDir alternation in `nextflow.config` (line 60) and drop
`BAM_TO_GFF` from it — the unfiltered GFF is now an intermediate. Full layout rework is
Phase C; for now `targets/<target>/filter_allmodels/` is acceptable.

## Checkpoint A

Run `nextflow run main.nf -profile test,crg --outdir results_15a`.

1. **Task counts** in `results_15a/pipeline_info/execution_trace.txt`: exactly **2**
   `MINIMAP2_ALIGN` (one per target), **2** `BAM_TO_GFF`, **2** `SPLIT_GFF_BY_SOURCE`,
   **1** `MERGE_SOURCE_FASTA`, **2** `TD2_PREDICT` in `PROJECTION` — down from 4 / 4 / — / — / 4.
2. **Per-source regression** — for each target × source:
   ```sh
   diff <(grep -v '^#' results_baseline/targets/T/bam_to_gff/T.from_S.lnc_RNA.projected.gff3) \
        <(grep -v '^#' results_15a/targets/T/.../T.from_S.lnc_RNA.projected.gff3)
   ```
   Expect **only deletions**, and every deleted block must be a whole gene/transcript/exon
   group whose transcript id appears in that run's TD2 `coding_ids`. A *modified* or *added*
   line means the merged alignment changed a per-read result — investigate, it should not
   happen (minimap2 is deterministic per read and `samtools sort` is stable, so within-source
   ordering is preserved).
3. **Byte-identity on the untouched tracks**: rerun with `--include_mrna --include_decoy` and
   confirm the `mRNA` and decoy per-source GFF3s are byte-identical to baseline.
4. **`allModels.raw.gff3` exists** per target × class, is grouped by seqid (`GFF_STATS` would
   exit 1 otherwise), and its transcript count equals the sum over the per-source splits.
5. **No duplicate-ID abort** on the test data; deliberately verify the guard fires by running
   the script by hand on a GFF3 with a duplicated `ID=`.
6. **Downstream unchanged**: `combined`, `top3`, `select_top_sources` and the gffcompare stats
   set have the same *filenames* as baseline, and `Lycaena_hippothoe_580924` still has no
   `gffcompare/`, no `select_top_sources/`. No `targets/null/` anywhere.
7. **Re-read against the plans**: confirm `GFFCOMPARE_COMBINE` and `consensus_top.nf` are
   still untouched (Phase A scope), that self-pairs are still excluded from the *old* combine
   (Phase B moves that), and that the meta contract in `CLAUDE.md` still holds —
   `feature_type`, `decoy`, `target_id`, `kind` present and `target_id` never null.

## Checkpoint A — RESULT (run 2026-08-11, `-profile test,crg`)

**PASSED.** Two runs: default (`results_15a`, 83 tasks) and `--include_mrna --include_decoy`
(`results_15a_full`, 235 tasks), each against an old-code baseline (`results_baseline`, and
`/nfs/scratch01/rg/fzanarello/fomo_results_baseline_full` produced from a `git worktree` at
`08e0d7a`).

| Check | Result |
|---|---|
| 1. Task counts, default | `MINIMAP2_ALIGN` 4→**2**, `BAM_TO_GFF` 4→**2**, `TD2_PREDICT` 4→**2**, `MERGE_SOURCE_FASTA` **1**, `SPLIT_GFF_BY_SOURCE` **2** |
| 1. Task counts, both flags | `MINIMAP2_ALIGN`/`BAM_TO_GFF`/`SAMTOOLS_STATS` 16→**8** (= F·D·T), `MERGE_SOURCE_FASTA` **4**, `TD2_PREDICT` 4→**2** |
| 2. Per-source regression | **all 4 byte-identical** |
| 3. Both flags | **all 16 byte-identical** — including `mRNA` and both decoy tracks |
| 4. `allModels.raw.gff3` | present per (target, class); transcripts = sum of splits (124, 125); seqids contiguous |
| 5. Duplicate-ID guard | verified by hand: exits 1 naming the id and both sources |
| 6. Downstream | `combined`, `top3`, every `.stats`, and `top_sources.csv` **content-identical** in both runs; no `targets/null`; the un-annotated target still has no `gffcompare/`, no `select_top_sources/`; 2 MultiQC reports |
| 7. Plan/contract re-read | `GFFCOMPARE_COMBINE` still N-input with self excluded; `consensus_top.nf` untouched; join keys intact; no nested subworkflows |

Three deviations from what this plan predicted, none of them failures:

- **Per-source files came out byte-identical, not "deletions only."** The pooled TD2 pass found
  **zero** coding transcripts on the test data (`n_coding=0` for both targets), so accepted
  consequence #1 in plan 15 did not materialise here. It stays structurally true — the split is
  downstream of the filter — but this data does not exercise it.
- **The gffcompare `.tmap`/`.refmap` filenames changed** (`…projected.filtered.gff3.tmap` →
  `…projected.gff3.tmap`), because gffcompare derives them from the query filename and the query
  is now the split output rather than a `.filtered.gff3`. `.stats`/`.tracking`/`.loci` are
  unchanged and `select_top_sources.py` parses only `.stats`, so nothing downstream breaks.
- **The projected TD2 pass is `T` tasks, not `D·T`** — it only ever runs on the real lncRNA
  class. Plan 15's cardinality table has been corrected.

Confirmed regression, as designed: the MultiQC samtools table collapses from one row per
(target, source) to one per (target, class) — `…allModels.1_lncRNA` with
`raw_total_sequences=145`, exactly the 64+81 of the two baseline rows, and `reads_mapped=136` =
64+72. `GFF_STATS_PROJECTED`'s per-source rows are unchanged.
