# Plan: FOMO Preprocessing — Decoy track (steps 02–04)

## Context

Step 00 is done: per source × feature_type (`lnc_RNA`, `mRNA`) the pipeline emits a filtered GFF3 and a spliced FASTA. The next preprocessing chunk produces **decoy** spliced FASTAs — source transcripts randomly relocated into the **target**'s intergenic intervals — which the projection stage uses as a negative-control track.

This plan covers legacy steps **02 → 03 → 04**. Step 01 (non-overlapping lncRNA filter) and step 05 (statistics) are explicitly deferred.

## Decisions

| Decision | Choice |
|----------|--------|
| Scope | Steps 02 (intergenic) + 03 (relocate) + 04 (decoy FASTA) |
| Step 02 tool | bedtools `complement` (drop AGAT) |
| Step 03 implementation | Refactor `legacy_scripts/preprocessing/03_relocate_loci.py` with argparse, then ship in `bin/` |
| Relocated feature types | Both `lnc_RNA` and `mRNA` |
| Source ↔ target pairing | Each source × single target (cartesian; current single-target model) |
| Decoy FASTA module | `include { GFFREAD as EXTRACT_DECOY_SEQUENCES }` (second alias) |
| Step 01 (non-overlap) | Skipped — revisit after projection results |
| Step 05 (statistics) | Out of scope for this PR |

## Channel topology

```
ch_input.target (singleton) [meta, fasta, gff3]
   │
   ├─ GUNZIP_TARGET_FASTA ─► [meta, fa]
   │       │
   │       ├─ SAMTOOLS_FAIDX ──► [meta, fai]
   │       │
   │       └────────────────────────────┐
   │                                    │
   └─ GFF_TO_GENE_BED ──► [meta, bed]   │
                  │                     │
                  └─ BEDTOOLS_COMPLEMENT(bed, fai) ──► ch_intergenic [meta_target, intergenic.bed]

FILTER_TRANSCRIPT.out.gff3 [meta(source,ftype), gff3]
                  ×
ch_intergenic (broadcast)
                  │
                  └─ RELOCATE_LOCI ──► [meta(source,ftype,decoy:true), decoy.gff3]
                                          │
                                          × GUNZIP_TARGET_FASTA.out.gunzip (broadcast)
                                          │
                                          └─ EXTRACT_DECOY_SEQUENCES ──► [meta, decoy.spliced.fasta]
```

Per run with N sources: 1 target gunzip + 1 faidx + 1 gene-BED + 1 complement + 2N relocate tasks + 2N decoy gffread tasks.

## Step 02 — Target intergenic intervals

Three modules, nf-core-style decomposition:

1. **`GUNZIP_TARGET_FASTA`** — reuse the existing `GUNZIP` nf-core module via a new alias (`include { GUNZIP as GUNZIP_TARGET_FASTA }`). Target FASTA needs decompressing once for both `SAMTOOLS_FAIDX` and `EXTRACT_DECOY_SEQUENCES`.
2. **`SAMTOOLS_FAIDX`** — nf-core module (`nf-core modules install samtools/faidx`). Output: `.fai` (chrom sizes derivable as `cut -f1,2`).
3. **`GFF_TO_GENE_BED`** — new local module (`modules/local/gff_to_gene_bed.nf`), ubuntu container. Two-line awk:
   ```bash
   gunzip -c ${gff3} 2>/dev/null || cat ${gff3} \
     | awk -F'\t' 'BEGIN{OFS="\t"} $0!~/^#/ && ($3=="gene" || $3=="pseudogene" || $3=="ncRNA_gene") { print $1, $4-1, $5 }' \
     | sort -k1,1 -k2,2n > ${prefix}.genes.bed
   ```
   No external deps beyond coreutils. Output: `[meta, bed]`.
4. **`BEDTOOLS_COMPLEMENT`** — nf-core module (`nf-core modules install bedtools/complement`). Inputs: sorted gene BED + genome file. Note: the nf-core module accepts a sizes file; the recipe in `modules.config` should set `ext.args2` (or pre-build a chrom-sizes file from `.fai`) — exact wiring to be confirmed against the installed module at implementation time.

Output: `ch_intergenic = [meta_target, intergenic.bed]`. Singleton, broadcast downstream via `.first()` or `combine`.

## Step 03 — Relocate source loci

### Script refactor (`bin/relocate_loci.py`)

Copied from legacy, with these changes:
- Replace positional argv with `argparse`:
  - `--input-gff` (required) — source filtered GFF3
  - `--intergenic-bed` (required) — target intergenic BED
  - `--output-gff` (required)
  - `--seed` (int, optional) — deterministic random placement
  - `--feature-type` (default `lnc_RNA`) — drives `transcript_type` in `read_gff3_models`
  - `--exon-type` (default `exon`)
  - `--min-intergenic-length` (int, default 0) — drop intervals shorter than this before random placement
- Keep all relocation logic (`relocate_model`, `remove_subinterval`, etc.) intact.
- Keep the `gene_decoy` / `*_decoy` ftype suffix on emitted features (legacy behavior — preserves downstream traceability).
- Print summary line to stderr (already present) and exit non-zero if zero models could be relocated (to surface bad inputs early).

### Module (`modules/local/relocate_loci.nf`)

```groovy
process RELOCATE_LOCI {
    tag "${meta.id}.${meta.feature_type}"
    label 'process_single'
    container 'python:3.11-slim'

    input:
    tuple val(meta), path(gff3), path(intergenic_bed)

    output:
    tuple val(meta), path("*.decoy.gff3"), emit: gff3

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.feature_type}.decoy"
    def seed   = params.relocate_seed ? "--seed ${params.relocate_seed}" : ''
    """
    relocate_loci.py \\
        --input-gff ${gff3} \\
        --intergenic-bed ${intergenic_bed} \\
        --output-gff ${prefix}.gff3 \\
        --feature-type ${meta.feature_type} \\
        ${seed} \\
        ${args}
    """
}
```

`params.relocate_seed = 42` in `nextflow.config` as the global default — overridable via `--relocate_seed` at runtime.

### Channel construction

```groovy
FILTER_TRANSCRIPT.out.gff3
    .combine( ch_intergenic.map { _meta, bed -> bed } )   // broadcast singleton
    .map { meta, gff3, bed -> tuple(meta + [decoy: true], gff3, bed) }
    .set { ch_relocate }

RELOCATE_LOCI(ch_relocate)
```

`meta` carries `[id, role:'source', feature_type, decoy:true]` from here on.

## Step 04 — Extract decoy spliced FASTA

Same shape as existing EXTRACT_SEQUENCES, with the alias `EXTRACT_DECOY_SEQUENCES` and the **target** fasta as the genome:

```groovy
include { GFFREAD as EXTRACT_DECOY_SEQUENCES } from '../../modules/nf-core/gffread/main'

RELOCATE_LOCI.out.gff3
    .combine( GUNZIP_TARGET_FASTA.out.gunzip.map { _m, fa -> fa } )  // broadcast
    .multiMap { meta, gff3, fa ->
        gff:   tuple(meta, gff3)
        fasta: fa
    }
    .set { ch_decoy_extract }

EXTRACT_DECOY_SEQUENCES(ch_decoy_extract.gff, ch_decoy_extract.fasta)
```

## `conf/modules.config` additions

```groovy
withName: 'FOMO:PREPROCESSING:GUNZIP_TARGET_FASTA' {
    ext.prefix = { "${meta.id}" }
}
withName: 'FOMO:PREPROCESSING:GFF_TO_GENE_BED' {
    ext.prefix = { "${meta.id}" }
}
withName: 'FOMO:PREPROCESSING:BEDTOOLS_COMPLEMENT' {
    ext.prefix = { "${meta.id}.intergenic" }
}
withName: 'FOMO:PREPROCESSING:RELOCATE_LOCI' {
    ext.prefix = { "${meta.id}.${meta.feature_type}.decoy" }
}
withName: 'FOMO:PREPROCESSING:EXTRACT_DECOY_SEQUENCES' {
    ext.args   = '-w'
    ext.prefix = { "${meta.id}.${meta.feature_type}.decoy.spliced" }
}
```

## Subworkflow integration

`subworkflows/local/preprocessing.nf` signature becomes:

```groovy
workflow PREPROCESSING {
    take:
    ch_sources   // [ meta, fasta, gff3 ]
    ch_target    // [ meta, fasta, gff3 ]  (singleton)

    main:
    // ... existing source path (GUNZIP_FASTA, FILTER_TRANSCRIPT, EXTRACT_SEQUENCES) ...

    // Target path → intergenic BED
    GUNZIP_TARGET_FASTA( ch_target.map { meta, fa, _g -> tuple(meta, fa) } )
    SAMTOOLS_FAIDX( GUNZIP_TARGET_FASTA.out.gunzip )
    GFF_TO_GENE_BED( ch_target.map { meta, _f, gff3 -> tuple(meta, gff3) } )
    BEDTOOLS_COMPLEMENT(
        GFF_TO_GENE_BED.out.bed,
        SAMTOOLS_FAIDX.out.fai.map { _m, fai -> fai }
    )

    // Decoy track
    // ... RELOCATE_LOCI + EXTRACT_DECOY_SEQUENCES per topology above ...

    emit:
    spliced_fasta       = EXTRACT_SEQUENCES.out.gffread_fasta
    filtered_gff3       = FILTER_TRANSCRIPT.out.gff3
    intergenic_bed      = BEDTOOLS_COMPLEMENT.out.bed
    decoy_gff3          = RELOCATE_LOCI.out.gff3
    decoy_spliced_fasta = EXTRACT_DECOY_SEQUENCES.out.gffread_fasta
}
```

`workflows/fomo.nf` already branches `ch_input.target` / `ch_input.source` — pass both into `PREPROCESSING`.

## Files to add

```
bin/relocate_loci.py                       # refactored from legacy 03_relocate_loci.py
modules/local/gff_to_gene_bed.nf
modules/local/relocate_loci.nf
modules/nf-core/samtools/faidx/main.nf     # via `nf-core modules install`
modules/nf-core/samtools/faidx/environment.yml
modules/nf-core/bedtools/complement/main.nf
modules/nf-core/bedtools/complement/environment.yml
```

Files to edit: `subworkflows/local/preprocessing.nf`, `workflows/fomo.nf`, `conf/modules.config`, `conf/test.config` (add `params.relocate_seed = 42`).

## Verification

1. `nextflow run . -profile test,docker` succeeds.
2. Outputs under `results/`:
   - 1× `<target>.intergenic.bed` (non-empty, sorted, no overlap with target gene intervals)
   - 8× `<source>.<ftype>.decoy.gff3` (4 sources × 2 ftypes)
   - 8× `<source>.<ftype>.decoy.spliced.fasta`
3. Sanity checks:
   - `bedtools intersect -a <target>.intergenic.bed -b <gene-bed> -u | wc -l` → 0
   - For each decoy GFF3: every transcript's `seqid` matches the target's chromosome names (not source).
   - Decoy spliced FASTA record count matches transcript count in matching decoy GFF3.
   - Re-running with the same `--seed` produces byte-identical decoy GFF3 (modulo header timestamps).
4. Compare a single decoy GFF3 against legacy output of `03_relocate_loci.py` on the same inputs+seed → structurally identical (exon block shapes preserved, only coordinates shifted).

## Deliberate omissions

- **Step 01** (non-overlapping lncRNA filter) — skipped, consistent with longest-isoform deferral.
- **Step 05** (AGAT statistics) — out of scope; will land in a follow-up alongside MultiQC integration.
- **Multiple targets** — single-target model retained; samplesheet still accepts exactly one `target` row.
- **Strand-aware relocation** — legacy script ignores strand collisions between decoy and existing features; we keep that behavior for parity.
