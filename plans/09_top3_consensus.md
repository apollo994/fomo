# Plan: Top-3 consensus annotation (selected by lncRNA transcript F1)

## Context

The pipeline already builds an **all-source** consensus per gene type
(`gffcompare` combine over every source projection) and benchmarks it against
the target reference. But pooling *every* source dilutes the consensus with
poorly-transferring annotations. This adds a second consensus that pools only
the **3 best sources**, where "best" is ranked by **lncRNA transcript-level F1**
(`F1 = 2·Sn·Pr/(Sn+Pr)`, from each source's gffcompare stats). Those same 3
sources feed both a lncRNA and an mRNA consensus. The result (`from_top3`) is
plugged into the same downstream as every other model — GTF→GFF3, AGAT,
gffcompare-vs-reference, MultiQC — so it can be compared head-to-head with the
per-source and all-source rows.

## Key constraint

Selection is **data-dependent on benchmark results**: the top-3 can only be
known after `BENCHMARKING` runs per-source gffcompare. So this is a **new
subworkflow after `BENCHMARKING`**, not an extension of `PROJECTION` (where the
all-source combine lives).

## Decisions (user-confirmed)

| Decision | Choice |
|----------|--------|
| Ranking metric | lncRNA **transcript-level** F1 only |
| What gets combined | lncRNA **and** mRNA projections from the same lncRNA-selected 3 sources |
| Decoys | none — real consensus only |
| Accuracy scatter | relocate `GFFCOMPARE_TO_SCATTER` to `REPORTING` so the `top3` dot shows over all stats |
| <3 sources | combine all available |
| F1 when Sn+Pr=0 | 0 (sinks empty/decoy-like projections) |

## Data flow

```
BENCHMARKING per-source lncRNA .stats ─► SELECT_TOP_SOURCES ─► top3.sources.csv (3 source IDs)
                                                                     │ splitCsv → [source]
PROJECTION.out.gff3, filtered to per-source REAL (!decoy && id!='combined')   (8 GFFs: 4 src × {lnc,mRNA})
   .map { meta,gff -> [meta.id, meta, gff] }
   .combine(top3_ids_keyed, by:0)            # inner join keeps the 3 selected sources → 6 GFFs
   .map  { _id, meta, gff -> [meta + [id:'top3'], gff, gtype] }
   .groupTuple(by gene type)                 # 2 groups: lnc_RNA, mRNA
        └─ GFFCOMPARE_COMBINE ─► 2 top3 consensus GTFs
             └─ COMBINED_GTF_TO_GFF (meta.id='top3') ─► GFF3
                  ├─ AGAT_SPSTATISTICS + AGAT_TO_MQC      (kind='projected', from_top3 rows)
                  └─ GFFCOMPARE vs target reference        (from_top3.<ft>.gffcompare.stats)
```

`from_top3` names fall out of the existing `ext.prefix` closures
(`${meta.target_id}.from_${meta.id}.${meta.feature_type}…`) because we set
`meta.id = 'top3'`, exactly like `from_combined` works today.

## New files
- `bin/gffcompare_stats.py` — small shared helper: parse a gffcompare `.stats`
  file → `{level: (sensitivity, precision)}` and derive `(source, ft, decoy)`
  from the filename. Refactored out of the duplicated logic in
  `bin/gffcompare_accuracy_mqc.py` (which then imports it).
- `bin/select_top_sources.py` — read all per-source lncRNA stats, compute
  transcript F1, emit the ranked top-N source IDs (default N=3) as a one-column
  CSV `source`.
- `modules/local/select_top_sources.nf` — wrapper (multiqc container, stdlib).
- `subworkflows/local/consensus_top.nf` — orchestration above.

## Edits to existing files
- `subworkflows/local/benchmarking.nf`
  - **emit** the filtered target references (`FILTER_TARGET.out.gff3`, keyed by
    feature_type) and keep emitting per-source `stats`.
  - **remove** the `GFFCOMPARE_TO_SCATTER` call + its `_mqc.json` from
    `mqc_files` (moves to REPORTING).
- `subworkflows/local/reporting.nf`
  - new `take:` input `ch_gffcompare_stats`; run `GFFCOMPARE_TO_SCATTER` over the
    union of benchmarking + consensus stats; add the json to the MultiQC inputs.
- `workflows/fomo.nf`
  - call `CONSENSUS_TOP(...)` after `BENCHMARKING`, passing per-source real
    projections, per-source lncRNA stats, target references, target fasta.
  - pass `BENCHMARKING.out.stats.mix(CONSENSUS_TOP.out.stats)` into `REPORTING`
    for the scatter; mix `CONSENSUS_TOP.out.mqc_files` into the report.
- `conf/modules.config`
  - `FOMO:CONSENSUS_TOP:GFFCOMPARE_COMBINE` → `ext.prefix = { "${meta.target_id}.top3.${meta.gtype}" }`
  - `FOMO:CONSENSUS_TOP:COMBINED_GTF_TO_GFF` → `from_top3.<ft>.combined`-style prefix
  - `FOMO:CONSENSUS_TOP:AGAT_SPSTATISTICS` / `AGAT_TO_MQC` → reuse projection prefix/section
  - `FOMO:CONSENSUS_TOP:GFFCOMPARE` → reuse benchmarking prefix/args (`-M --no-merge`)
  - `FOMO:CONSENSUS_TOP:SELECT_TOP_SOURCES` → (optional) `ext.args` for N via `params.top_consensus_n`
- `nextflow.config` + `nextflow_schema.json`
  - new param `top_consensus_n = 3` (kept in sync per CLAUDE.md).

## `select_top_sources.py` sketch
```
in : N lncRNA .stats files + --n 3
for each: source = parse_identity(); sn,pr = parse_accuracy()['Transcript']
          f1 = 0 if sn+pr==0 else 2*sn*pr/(sn+pr)
rank sources by f1 desc; emit top-N source ids → CSV (header 'source')
```

## Selection → filter mechanics (Nextflow)
```groovy
ch_top = SELECT_TOP_SOURCES.out.csv
    .splitCsv(header: true)
    .map { row -> tuple(row.source, true) }          // [source, _]

ch_persource_real = ch_projected_real                // !decoy && id != 'combined'
    .map { meta, gff -> tuple(meta.id, meta, gff) }

ch_selected = ch_persource_real
    .combine(ch_top, by: 0)                          // inner join on source id
    .map { _src, meta, gff, _ ->
        def gtype = meta.feature_type                // lnc_RNA | mRNA (no decoys here)
        tuple([id: 'top3', target_id: meta.target_id, feature_type: gtype,
               decoy: false, gtype: gtype], gff)
    }
    .groupTuple()
    .map { meta, gffs -> tuple(meta, gffs.sort { it.name }) }
GFFCOMPARE_COMBINE(ch_selected, [[:],[],[]], [[:],[]])
```

## Verification
1. `nextflow run main.nf -profile test,docker` completes.
2. `SELECT_TOP_SOURCES` emits exactly 3 source IDs; cross-check by hand against
   the lncRNA transcript Sn/Pr in `results/gffcompare/*from_*.lnc_RNA.*.stats`
   (top-3 F1).
3. `results/.../*.top3.lnc_RNA.combined.gtf` and `*.top3.mRNA.combined.gtf` exist;
   each combined from exactly the 3 selected sources (check `num_samples` /
   `.tracking` column count = 3).
4. AGAT projection table + benchmarking gffcompare both show `from_top3.lnc_RNA`
   and `from_top3.mRNA` rows.
5. Custom accuracy scatter includes a `top3` dot (distinct colour) for lncRNA and
   mRNA, alongside per-source and `combined`.
6. `-resume` is stable.

## Out of scope
- top-3 decoy consensus (user: real only).
- mRNA-based or global ranking (user: lncRNA transcript F1 only).
- Tunable F1 level (fixed to transcript).
