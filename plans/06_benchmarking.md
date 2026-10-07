# Plan: FOMO Benchmarking — gffcompare per projection

## Context

From `BRAINSTORM.md` benchmarking step:

```
- run gffcompare (gffcompare -M --no-exon-merge)
- extract gffcompare statistics for report
```

The pipeline currently produces, per target species, **16 projected GFF3s** (4 sources × 4 classes: `mRNA`, `lnc_RNA`, `mRNA.decoy`, `lnc_RNA.decoy`). Benchmarking compares each projection against the target species' real annotation so we can quantify sensitivity/precision per source and per feature class — and use the decoy projections as a false-positive baseline.

No legacy reference script — `BRAINSTORM.md` directly prescribes the gffcompare invocation.

## Decisions

| Decision | Choice |
|----------|--------|
| Reference scope | **Per-class filtered target.** `mRNA`/`mRNA.decoy` ↔ target mRNA; `lnc_RNA`/`lnc_RNA.decoy` ↔ target lnc_RNA. Apples-to-apples. |
| Decoy benchmarking | **Option A — same-class real annotation.** Free byproduct of per-class filtering: decoys pair on `feature_type` like real projections. Sp on decoy rows IS the false-positive baseline; "true precision lift" = non-decoy Sp − decoy Sp. |
| `-M` (discard single-exon transfrags) | **Run twice.** Two parallel `GFFCOMPARE` aliases: `GFFCOMPARE_M` with `-M --no-exon-merge`, `GFFCOMPARE_NOM` with `--no-exon-merge`. Lets us see both the multi-exon view (BRAINSTORM-prescribed) and the full view (where single-exon decoys/lncRNAs are visible). 32 tasks total. |
| Genome FASTA into gffcompare | **Omit** — optional input, no sequence-level metrics needed. |
| Reuse `FILTER_TRANSCRIPT` | **Yes** — same module already used for sources. Aliased as `FILTER_TARGET` so target-specific `ext.prefix` doesn't collide with the source-side selector. |
| Output integration | Pipe `.stats` files into the existing **MULTIQC** call in `STATS`. MultiQC has a native gffcompare module — no custom TSV parser needed. |
| Subworkflow | New `subworkflows/local/benchmarking.nf`. Owns target filtering + both GFFCOMPARE aliases. |
| Output naming | `<target>.from_<source>.<ftype>[.decoy].M.gffcompare.*` and `<target>.from_<source>.<ftype>[.decoy].noM.gffcompare.*` |

## 1. Install nf-core/gffcompare

```
nf-core modules install gffcompare
```

(If the CLI hits the same TTY issue as before, fetch `main.nf` + `environment.yml` from the master branch manually — same workflow used for `minimap2/index` and `minimap2/align`.)

Module signature:

```
input:
  tuple val(meta),  path(gtfs)             // query GFF/GTF(s)
  tuple val(meta2), path(fasta), path(fai) // optional — pass [[], [], []]
  tuple val(meta3), path(reference_gtf)    // reference annotation
output:
  tuple val(meta), path("*.stats"),         emit: stats
  tuple val(meta), path("*.annotated.gtf"), emit: annotated_gtf, optional
  tuple val(meta), path("*.loci"),          emit: loci
  tuple val(meta), path("*.tracking"),      emit: tracking
  tuple val(meta), path("*.tmap"),          emit: tmap, optional
  tuple val(meta), path("*.refmap"),        emit: refmap, optional
```

## 2. `BENCHMARKING` subworkflow

`subworkflows/local/benchmarking.nf`:

```groovy
include { GUNZIP as GUNZIP_TARGET_GFF             } from '../../modules/nf-core/gunzip/main'
include { FILTER_TRANSCRIPT as FILTER_TARGET      } from '../../modules/local/filter_transcript'
include { GFFCOMPARE as GFFCOMPARE_M              } from '../../modules/nf-core/gffcompare/main'
include { GFFCOMPARE as GFFCOMPARE_NOM            } from '../../modules/nf-core/gffcompare/main'

workflow BENCHMARKING {
    take:
    ch_target          // [ meta(role:'target'), fasta.gz, gff3.gz ] × 1
    ch_projected_gff3  // [ meta(target_id, id, feature_type, decoy), gff3 ] × 16

    main:

    GUNZIP_TARGET_GFF(
        ch_target.map { meta, _fa, gff3 -> tuple(meta, gff3) }
    )

    // Filter target by feature_type — produces target.mRNA.gff3 + target.lnc_RNA.gff3
    ch_target_filter_in = GUNZIP_TARGET_GFF.out.gunzip
        .combine(Channel.of('lnc_RNA', 'mRNA'))
        .map { meta, gff3, ftype ->
            tuple(meta + [feature_type: ftype, decoy: false], gff3, ftype)
        }

    FILTER_TARGET(ch_target_filter_in)

    // Pair each projection with the target-class reference of the same feature_type
    ch_paired = ch_projected_gff3
        .map { meta, gff3 -> tuple(meta.feature_type, meta, gff3) }
        .combine(
            FILTER_TARGET.out.gff3.map { meta, gff3 -> tuple(meta.feature_type, gff3) },
            by: 0
        )

    ch_query = ch_paired.multiMap { _ftype, q_meta, q_gff, _r_gff ->
        with_meta: tuple(q_meta, q_gff)
    }

    ch_empty_ref = ch_paired.multiMap { _ftype, q_meta, _q_gff, _r_gff ->
        empty: tuple([id: q_meta.target_id], [], [])
    }

    ch_reference = ch_paired.multiMap { _ftype, q_meta, _q_gff, r_gff ->
        ref: tuple([id: "${q_meta.target_id}.${q_meta.feature_type}"], r_gff)
    }

    GFFCOMPARE_M(
        ch_query.with_meta,
        ch_empty_ref.empty,
        ch_reference.ref
    )

    GFFCOMPARE_NOM(
        ch_query.with_meta,
        ch_empty_ref.empty,
        ch_reference.ref
    )

    emit:
    stats_m   = GFFCOMPARE_M.out.stats       // [ meta, *.stats ] × 16
    stats_nom = GFFCOMPARE_NOM.out.stats     // [ meta, *.stats ] × 16
}
```

> **Channel-topology gotcha:** `ch_paired` is a queue channel; using it three times via `multiMap` keeps the 16-element ordering intact for all three branches and the cross-branch pairing required by `GFFCOMPARE`. Don't `.first()` — the reference per task is class-specific, not a single broadcast value.

## 3. Wire into top-level

`workflows/fomo.nf`:

```groovy
include { BENCHMARKING } from '../subworkflows/local/benchmarking'

...
PROJECTION(
    ch_input.target,
    PREPROCESSING.out.spliced_fasta,
    PREPROCESSING.out.decoy_spliced_fasta
)

BENCHMARKING(ch_input.target, PROJECTION.out.gff3)

STATS(
    ch_input.source,
    ch_input.target,
    PREPROCESSING.out.filtered_gff3,
    PREPROCESSING.out.decoy_gff3,
    PREPROCESSING.out.spliced_fasta,
    PREPROCESSING.out.decoy_spliced_fasta,
    BENCHMARKING.out.stats_m
        .mix(BENCHMARKING.out.stats_nom)        // mix both M modes into the MultiQC bundle
)
```

## 4. `STATS` subworkflow — accept gffcompare stats

`subworkflows/local/stats.nf` — add the new `take:` input and mix it into `ch_mqc_files`:

```groovy
take:
    ...
    ch_gffcompare_stats  // [ meta, *.stats ] × 32  (16 with-M + 16 no-M)

main:
    ...
    ch_mqc_files = SEQKIT_TO_MQC.out.mqc.map { _m, t -> t }
        .mix(AGAT_TO_MQC.out.summary.map { _m, t -> t })
        .mix(AGAT_TO_MQC.out.full   .map { _m, t -> t })
        .mix(ch_gffcompare_stats     .map { _m, s -> s })   // NEW
        .collect()
```

MultiQC's native `gffcompare` module auto-detects `*.stats`. No custom-content config addition needed.

## 5. `assets/multiqc_config.yml` tweak

Add `.gffcompare` and the `.M`/`.noM` suffix to `extra_fn_clean_exts` so MultiQC sample labels collapse `<target>.from_<source>.<ftype>[.decoy].M.gffcompare.stats` → `<target>.from_<source>.<ftype>[.decoy].M`:

```yaml
extra_fn_clean_exts:
  ...
  - ".gffcompare"
```

(`.stats` is already in the list.) The `.M` / `.noM` distinguishes the two gffcompare runs in the MultiQC table — that's intentional, not stripped.

## 6. `conf/modules.config` additions

```groovy
withName: 'FOMO:BENCHMARKING:GUNZIP_TARGET_GFF' {
    ext.prefix = { "${meta.id}.raw" }
}

withName: 'FOMO:BENCHMARKING:FILTER_TARGET' {
    ext.prefix = { "${meta.id}.${meta.feature_type}" }
}

withName: 'FOMO:BENCHMARKING:GFFCOMPARE_M' {
    ext.args   = '-M --no-exon-merge'
    ext.prefix = { "${meta.target_id}.from_${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.M.gffcompare" }
}

withName: 'FOMO:BENCHMARKING:GFFCOMPARE_NOM' {
    ext.args   = '--no-exon-merge'
    ext.prefix = { "${meta.target_id}.from_${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.noM.gffcompare" }
}
```

## 7. Expected outputs

**32 gffcompare runs** (16 projections × 2 M modes). Each produces:

```
<target>.from_<source>.<ftype>[.decoy].{M,noM}.gffcompare.stats
<target>.from_<source>.<ftype>[.decoy].{M,noM}.gffcompare.annotated.gtf
<target>.from_<source>.<ftype>[.decoy].{M,noM}.gffcompare.loci
<target>.from_<source>.<ftype>[.decoy].{M,noM}.gffcompare.tracking
<target>.from_<source>.<ftype>[.decoy].{M,noM}.gffcompare.<query>.tmap
<target>.from_<source>.<ftype>[.decoy].{M,noM}.gffcompare.<query>.refmap
```

For the current 4-source test samplesheet, 32 rows appear in the MultiQC `gffcompare` section — 16 with-M, 16 without-M, each labeled by source, feature class, decoy/not.

## 8. Verification

1. `nextflow run . -profile test,docker -resume`:
   - 1 GUNZIP_TARGET_GFF + 2 FILTER_TARGET + 16 GFFCOMPARE_M + 16 GFFCOMPARE_NOM = **35 new tasks**.
   - STATS re-runs because its `take:` signature changes (forces MULTIQC re-run too — expected).
2. Sanity checks (open the MultiQC report):
   - **Real projections, with-M:** mRNA Sn highest from closely-related sources; mRNA Sp should be moderate-to-high; lncRNA Sn lower (small N) but >0; lncRNA Sp moderate.
   - **Decoy projections, with-M:** Sn ≈ 0; Sp gives same-class FPR — should be much lower than non-decoy Sp.
   - **No-M runs:** generally higher counts because single-exon transfrags re-enter; biggest impact on decoys and lncRNAs (most often single-exon).
3. "True precision lift" sanity: for at least one (source, feature_type) pair, `Sp(real) − Sp(decoy)` should be a clearly positive number.
4. File hygiene: `ls results/gffcompare* | wc -l` matches 32× artifact count (stats + annotated.gtf + loci + tracking + tmap + refmap ≈ 192 files).

## 9. Deliberate omissions

- **Cross-class benchmarking** (mRNA proj vs target lncRNA, etc.): out of scope; would dilute the class-vs-class signal.
- **Whole-target reference for decoys**: rejected — same-class is the cleaner FPR baseline; per the user's confirmation.
- **Aggregating Sn/Sp across sources** into a single per-target summary: useful but a separate plan; MultiQC table already supports per-row sorting.
- **gffcompare `--combined-gtf` multi-query mode**: not used; we want per-projection signal.
- **Sequence-level (`-r ref.fa`) metrics**: not needed for these accuracy stats.
