# Plan: FOMO Mapping — Project spliced FASTAs onto target

## Context

Preprocessing (subworkflow `PREPROCESSING`) currently produces, per source species, **three** spliced FASTAs:

- `mRNA` spliced exons (`meta.decoy == false`, `meta.feature_type == 'mRNA'`)
- `lnc_RNA` spliced exons (`meta.decoy == false`, `meta.feature_type == 'lnc_RNA'`)
- decoy spliced exons (`meta.decoy == true`, `meta.feature_type ∈ {mRNA, lnc_RNA}`)

The next stage from `BRAINSTORM.md` is **projection** — splice-aware alignment of every source spliced FASTA onto the single target genome with minimap2, producing one sorted+indexed BAM per (source, feature_type, decoy?) combination. Reference implementation: `legacy_scripts/minimap_transfer/00_run_minimap_base.sh`.

This plan covers **only the FASTA → BAM mapping step** (legacy script `00`). Conversion BAM → GFF3 (legacy `01_convert_bam_to_gff.sh`) and metric extraction/plots (legacy `02`–`04`) are explicitly out of scope and will come in follow-up PRs.

## Decisions

| Decision | Choice |
|----------|--------|
| Alignment tool | `minimap2 -ax splice -uf` (`-uf` forces forward-strand reads — appropriate for spliced transcript exons; deviates from legacy `00_run_minimap_base.sh` which omitted it) |
| Reference handling | Build a `.mmi` index **once** with `nf-core/minimap2/index`, then reuse across every alignment. Avoids re-indexing the same target FASTA per (source × feature_type × decoy) task. |
| nf-core modules | `nf-core/minimap2/index` + `nf-core/minimap2/align` (align wraps `minimap2 \| samtools sort \| samtools index`, emits sorted BAM + BAI) |
| BAM index | Yes (`bam_index_extension = 'bai'`) — needed downstream by validation/projection conversion |
| Output format | Sorted BAM (legacy script's intermediate SAM is discarded inside the module) |
| samtools stats | **Skipped here** — produced later in a dedicated stats step alongside the existing MultiQC report |
| Channel cardinality | One BAM per (source × feature_type × decoy?) — for current 4-source samplesheet that's ~16 BAMs (4 × {mRNA, lncRNA, mRNA-decoy, lncRNA-decoy}); confirm against `RENAME_FASTA_HEADERS` outputs at impl time |

## 1. New nf-core modules to install

```
nf-core modules install minimap2/index
nf-core modules install minimap2/align
```

This drops `modules/nf-core/minimap2/{index,align}/` (each with `main.nf` + `meta.yml` + `tests/`). No local module needed for this step.

Module signatures (from `describe_nfcore_module`):

```
# minimap2/index
input:  tuple val(meta), path(fasta)
output: tuple val(meta), path("*.mmi"), emit: index

# minimap2/align
input:
  tuple val(meta),  path(reads)        // the source spliced FASTA
  tuple val(meta2), path(reference)    // .mmi index from minimap2/index (FASTA also accepted)
  val   bam_format                     // true  → emit BAM (not PAF)
  val   bam_index_extension            // 'bai' → also emit *.bam.bai
  val   cigar_paf_format               // false
  val   cigar_bam                      // false
```

## 2. Target FASTA preparation

Target compressed FASTA enters via `ch_input.target` in `workflows/fomo.nf`. Steps:

1. **Gunzip** the target FASTA once with `GUNZIP` (already vendored under `modules/nf-core/gunzip`) — matches the source-side pattern in `PREPROCESSING`.
2. **Index** the decompressed target FASTA once with `MINIMAP2_INDEX`, passing `ext.args = '-x splice'` so the `.mmi` is built with the same preset used at alignment time (otherwise minimap2 warns / silently changes parameters when align preset differs from index preset).
3. **Broadcast** the resulting `.mmi` as a value channel into every `MINIMAP2_ALIGN` call.

## 3. New subworkflow: `MAPPING`

`subworkflows/local/mapping.nf`:

```groovy
include { GUNZIP as GUNZIP_TARGET } from '../../modules/nf-core/gunzip/main'
include { MINIMAP2_INDEX          } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN          } from '../../modules/nf-core/minimap2/align/main'

workflow MAPPING {
    take:
    ch_target            // [ meta(role:'target'), fasta.gz, gff3.gz ] × 1
    ch_spliced_fasta     // [ meta(decoy:false, feature_type), fasta ] × 2N  (from PREPROCESSING)
    ch_decoy_spliced     // [ meta(decoy:true,  feature_type), fasta ] × 2N  (from PREPROCESSING)

    main:

    // 1. Decompress the target FASTA once
    GUNZIP_TARGET( ch_target.map { meta, fa, _gff3 -> tuple(meta, fa) } )

    // 2. Build .mmi index once (uses ext.args = '-x splice' from modules.config)
    MINIMAP2_INDEX( GUNZIP_TARGET.out.gunzip )

    // 3. Broadcast index as value channel
    ch_reference = MINIMAP2_INDEX.out.index
        .map { meta, mmi -> tuple([id: meta.id], mmi) }
        .first()

    // 4. Reads channel: every source spliced FASTA (mRNA / lncRNA / decoy).
    //    Attach target id to meta so output filenames encode both target and source.
    ch_reads = ch_spliced_fasta
        .mix(ch_decoy_spliced)
        .combine(ch_reference)
        .map { meta, reads, ref_meta, _mmi ->
            tuple(meta + [target_id: ref_meta.id], reads)
        }

    MINIMAP2_ALIGN(
        ch_reads,
        ch_reference,
        true,      // bam_format
        'bai',     // bam_index_extension
        false,     // cigar_paf_format
        false      // cigar_bam
    )

    emit:
    bam      = MINIMAP2_ALIGN.out.bam       // [ meta, *.bam ]
    bai      = MINIMAP2_ALIGN.out.index     // [ meta, *.bam.bai ]
    versions = MINIMAP2_INDEX.out.versions.mix(MINIMAP2_ALIGN.out.versions)
}
```

> **meta convention:** the read-side `meta` carries `id` (source species), `feature_type`, and `decoy`. We add `target_id` so per-output filenames include the target — e.g. `Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.lnc_RNA.spliced.bam`. This flattens the legacy directory layout `<target>_from_<source>/{pcg,lnc,decoy}/aligned.bam` into a single filename.

## 4. `conf/modules.config` additions

```groovy
withName: 'FOMO:MAPPING:GUNZIP_TARGET' {
    ext.prefix = { "${meta.id}" }
}

withName: 'FOMO:MAPPING:MINIMAP2_INDEX' {
    ext.args   = '-x splice'   // MUST match the preset used at align time, otherwise minimap2 warns / silently overrides
    ext.prefix = { "${meta.id}.splice" }
    cpus   = 2
    memory = 8.GB
    time   = 30.min
}

withName: 'FOMO:MAPPING:MINIMAP2_ALIGN' {
    // Module template prepends `-a` when bam_format=true → reproduces `minimap2 -ax splice -uf`.
    // -uf forces forward-strand reads (correct for spliced transcript exons).
    ext.args  = '-x splice -uf'
    ext.prefix = { "${meta.target_id}.from_${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.spliced" }
    cpus   = 4
    memory = 8.GB
    time   = 2.h
}
```

> Verify the `-a` injection against the installed module template at implementation time; if the template doesn't inject `-a`, change `ext.args` to `'-ax splice -uf'`.

## 5. Wire into `workflows/fomo.nf`

```groovy
include { MAPPING } from '../subworkflows/local/mapping'

workflow FOMO {
    ...
    PREPROCESSING(ch_input.source)
    STATS(...)

    MAPPING(
        ch_input.target,
        PREPROCESSING.out.spliced_fasta,
        PREPROCESSING.out.decoy_spliced_fasta
    )
}
```

## 6. Expected outputs

For the current 4-source / 1-target test samplesheet:

- **16 BAMs** = 4 sources × (mRNA + lncRNA + mRNA-decoy + lncRNA-decoy)
- **16 BAIs** alongside

Files land under `results/minimap2/align/` (default nf-core publishDir) with names like:

```
Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.mRNA.spliced.bam
Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.mRNA.spliced.bam.bai
Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.lnc_RNA.spliced.bam
Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.mRNA.decoy.spliced.bam
Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.lnc_RNA.decoy.spliced.bam
...
```

(Cardinality assumes one decoy FASTA per (source, feature_type) as currently emitted by `PREPROCESSING.out.decoy_spliced_fasta` — confirm against `RENAME_FASTA_HEADERS` outputs in case decoy is emitted only per source rather than per feature_type.)

## 7. Verification

1. `nextflow run . -profile test,docker -resume` — confirm `MAPPING` runs after `PREPROCESSING` without retriggering it.
2. Check that the expected number of sorted BAMs + BAIs are produced.
3. Spot-check one BAM against the legacy script on the same input:

   ```
   samtools view -c -F 2308 results/...mRNA.spliced.bam
   ```

   Compare to the count from running `legacy_scripts/minimap_transfer/00_run_minimap_base.sh` directly on the same `<target.fa, source.spliced.fa>` pair.
4. `samtools quickcheck` on every BAM.

## 8. Deliberate omissions (follow-up PRs)

- **BAM → GFF3 conversion** (legacy `01_convert_bam_to_gff.sh`). Needs a new local module wrapping the two-pass `samtools view | bedtools bamtobed -bed12 | awk` pipeline that injects `NM/AS/de` tags onto transcript rows. Track as `plans/04_bam_to_gff.md`.
- **Alignment QC** — samtools stats + the per-feature-class summaries the legacy script writes (`*.stats.all`, `*.stats.primary`). Will be folded into the existing MultiQC report in a stats follow-up.
- **`02_extract_alignment_metrics.py` / `03_plot_alignment_metrics.py`** — NM/AS/de aggregation + plots. Out of scope until BAM→GFF3 lands.
- **Longest-isoform / non-overlapping lncRNA filtering** — still skipped (per `plans/00_preprocessing-step00.md`); revisit when reviewing projection quality.
