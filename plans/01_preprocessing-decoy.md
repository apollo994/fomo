# Plan: FOMO Preprocessing — Decoy track (steps 02–04)

## Context

Step 00 is done: per source × feature_type (`lnc_RNA`, `mRNA`) the pipeline emits a filtered GFF3 and a spliced FASTA. The next preprocessing chunk produces **decoy** spliced FASTAs — source transcripts randomly relocated into **each source's own** intergenic intervals — which the projection stage uses as a negative-control track. The decoy sequences are extracted from the source FASTA (not the target) because only source-derived sequences are aligned to the target during projection; decoys must mimic source lncRNA composition to serve as a meaningful FP control.

This plan covers legacy steps **02 → 03 → 04**. Step 01 (non-overlapping lncRNA filter) and step 05 (statistics) are explicitly deferred.

## Decisions

| Decision | Choice |
|----------|--------|
| Scope | Steps 02 (intergenic) + 03 (relocate) + 04 (decoy FASTA) |
| Step 02 tool | bedtools `complement` on a sorted+merged gene BED (drop AGAT) |
| Step 02 placement | Per **source** species (N runs), not target — decoys serve as source-side FP controls |
| Step 03 implementation | Refactor `legacy_scripts/preprocessing/03_relocate_loci.py` with argparse, then ship in `bin/` |
| Relocated feature types | Both `lnc_RNA` and `mRNA` |
| Chromosome selection | Weighted by total eligible intergenic length per chromosome (not uniform) |
| Decoy marker | `decoy=true` GFF3 attribute on gene + transcript (NOT a `*_decoy` ftype suffix — gffread doesn't recognize non-standard transcript types) |
| Empty-input behavior | Emit empty GFF3 with header + warning to stderr (exit 0) — pipeline keeps running when test data has no models on a chromosome |
| Decoy FASTA module | `include { GFFREAD as EXTRACT_DECOY_SEQUENCES }` (second alias); uses source FASTA |
| Step 01 (non-overlap) | Skipped — revisit after projection results |
| Step 05 (statistics) | Out of scope for this PR |

## Channel topology

```
ch_sources [meta, fasta, gff3] × N
   │
   ├─ GUNZIP_FASTA ──► [meta(source), fa]
   │       │
   │       └─ SAMTOOLS_FAIDX(fa, get_sizes=true) ──► [meta, sizes]
   │
   └─ GFF_TO_GENE_BED ──► [meta(source), sorted+merged genes.bed]
                              │
                              × (join by meta.id) sizes
                              │
                              └─ BEDTOOLS_COMPLEMENT ──► ch_intergenic [meta(source), intergenic.bed]

FILTER_TRANSCRIPT.out.gff3 [meta(source, ftype, decoy:false), gff3]
                  × (join by meta.id) ch_intergenic
                  │
                  └─ RELOCATE_LOCI ──► [meta(source, ftype, decoy:true), decoy.gff3]
                                              │
                                              × (join by meta.id) GUNZIP_FASTA.out.gunzip
                                              │
                                              └─ EXTRACT_DECOY_SEQUENCES ──► [meta, decoy.spliced.fasta]

EXTRACT_SEQUENCES.out.gffread_fasta ──┐
                                       ├─ mix ─► RENAME_FASTA_HEADERS ──► [meta, *.renamed.fasta]
EXTRACT_DECOY_SEQUENCES.out ──────────┘                                   (split by meta.decoy on emit)
```

Per run with N sources: N gunzip + N faidx + N gene-BED + N complement + 2N filter + 2N extract + 2N relocate + 2N decoy gffread + 4N rename.

## Step 02 — Per-source intergenic intervals

Three modules, nf-core-style decomposition, **run once per source species**:

1. **`SAMTOOLS_FAIDX`** — nf-core module. Indexes each decompressed source FASTA (reusing `GUNZIP_FASTA.out.gunzip` from step 00) with `get_sizes = true` to emit a `.sizes` file (`cut -f 1,2` of `.fai`).
2. **`GFF_TO_GENE_BED`** — local module (`modules/local/gff_to_gene_bed.nf`), ubuntu container. Pipeline runs with `set -euo pipefail`:
   ```bash
   gunzip -c "${gff3}" | awk '... gene/pseudogene/ncRNA_gene → BED3 ...' \
       | sort -k1,1 -k2,2n \
       | awk '... merge overlapping/adjacent intervals ...' \
       > ${prefix}.genes.bed
   ```
   The inline awk-based merge replaces a separate `bedtools merge` call; required because Ensembl annotations have overlapping gene/pseudogene/ncRNA_gene records that would otherwise break `bedtools complement`.
3. **`BEDTOOLS_COMPLEMENT`** — nf-core module. Inputs: sorted+merged gene BED + `.sizes` file. Output: per-source intergenic BED.

Channel construction pairs `GFF_TO_GENE_BED.out.bed` with `SAMTOOLS_FAIDX.out.sizes` by `meta.id` using `combine(by:0)` + `multiMap` to keep emission pairing safe.

## Step 03 — Relocate source loci

### Script refactor (`bin/relocate_loci.py`)

Copied from legacy, with these changes:
- Replace positional argv with `argparse`:
  - `--input-gff` (required) — source filtered GFF3
  - `--intergenic-bed` (required) — source-species intergenic BED
  - `--output-gff` (required)
  - `--seed` (int, optional) — deterministic random placement
  - `--feature-type` (default `lnc_RNA`) — drives `transcript_type` in `read_gff3_models`
  - `--exon-type` (default `exon`)
  - `--min-intergenic-length` (int, default 0) — drop intervals shorter than this before random placement
- Keep all relocation logic (`relocate_model`, `remove_subinterval`, etc.) intact.
- Chromosome placement: legacy script required `source_chrom == target_chrom`, which always fails cross-species. Replaced with weighted random selection across all source chromosomes that have at least one interval ≥ transcript span (weight = total eligible intergenic length on that chromosome, so big chromosomes receive proportionally more decoys).
- **Decoy marker**: emit standard transcript ftypes (`lnc_RNA`, `mRNA`) and standard `gene` ftype, with a `decoy=true` attribute on both gene and transcript records. (The original plan called for `*_decoy` ftypes, but gffread does not recognize non-standard transcript types and silently emits zero records.)
- Print summary line to stderr (already present). On zero models relocated: emit a valid empty GFF3 (header only) and warn to stderr; do NOT exit non-zero, since with subsampled test data a single chromosome may legitimately contain zero transcripts of a given feature type.

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

Same shape as existing EXTRACT_SEQUENCES, with the alias `EXTRACT_DECOY_SEQUENCES` and the **source** fasta as the genome (the decoy GFF3's seqids are source chromosome names; sequence content must come from the source assembly):

```groovy
include { GFFREAD as EXTRACT_DECOY_SEQUENCES } from '../../modules/nf-core/gffread/main'

RELOCATE_LOCI.out.gff3
    .map { meta, gff3 -> tuple(meta.id, meta, gff3) }
    .combine( GUNZIP_FASTA.out.gunzip.map { m, fa -> tuple(m.id, fa) }, by: 0 )
    .multiMap { id, meta, gff3, fa ->
        gff:   tuple(meta, gff3)
        fasta: fa
    }
    .set { ch_decoy_extract }

EXTRACT_DECOY_SEQUENCES(ch_decoy_extract.gff, ch_decoy_extract.fasta)
```

## Step 05 — Rename FASTA headers

After both extract steps, a local `RENAME_FASTA_HEADERS` module rewrites each FASTA header to a pipe-delimited form `<transcript_id>|<type>|<species>` so downstream alignment hits can be traced back to source transcript + type + species. The same module handles source and decoy FASTAs — `meta.decoy` and `meta.feature_type` drive the type label (`lncRNA`, `mRNA`, `decoy_lncRNA`, `decoy_mRNA`). Mix the two extract outputs into one channel and split back on `emit` via `.filter { meta, _ -> meta.decoy }`.

## `conf/modules.config` additions

```groovy
withName: 'FOMO:PREPROCESSING:GFF_TO_GENE_BED' {
    ext.prefix = { "${meta.id}" }
}
withName: 'FOMO:PREPROCESSING:BEDTOOLS_COMPLEMENT' {
    ext.prefix = { "${meta.id}.intergenic" }
}
withName: 'FOMO:PREPROCESSING:RELOCATE_LOCI' {
    ext.prefix = { "${meta.id}.${meta.feature_type}" }
}
withName: 'FOMO:PREPROCESSING:EXTRACT_DECOY_SEQUENCES' {
    ext.args   = '-w'
    ext.prefix = { "${meta.id}.${meta.feature_type}.decoy.spliced" }
}
withName: 'FOMO:PREPROCESSING:RENAME_FASTA_HEADERS' {
    ext.prefix = { "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.spliced" }
}
```

## Subworkflow integration

`subworkflows/local/preprocessing.nf` signature stays as `take: ch_sources` — the target is not needed for the decoy track (target FASTA is consumed in later projection stages). Per-source intergenic + decoy + rename steps slot in after the existing source path.

Emit channels:

```groovy
emit:
spliced_fasta       = RENAME_FASTA_HEADERS.out.fasta.filter { meta, fa -> !meta.decoy }
filtered_gff3       = FILTER_TRANSCRIPT.out.gff3
intergenic_bed      = BEDTOOLS_COMPLEMENT.out.bed
decoy_gff3          = RELOCATE_LOCI.out.gff3
decoy_spliced_fasta = RENAME_FASTA_HEADERS.out.fasta.filter { meta, fa -> meta.decoy }
```

The `decoy: false` flag is set explicitly on the source path (in the same `.map` that adds `feature_type`) so downstream filters can rely on `meta.decoy` being a proper boolean rather than relying on null-truthiness.

## Files to add

```
bin/relocate_loci.py                       # refactored from legacy 03_relocate_loci.py
modules/local/gff_to_gene_bed.nf
modules/local/relocate_loci.nf
modules/local/rename_fasta_headers.nf
modules/nf-core/samtools/faidx/main.nf     # fetched from nf-core/modules master
modules/nf-core/samtools/faidx/environment.yml
modules/nf-core/bedtools/complement/main.nf
modules/nf-core/bedtools/complement/environment.yml
```

Files to edit: `subworkflows/local/preprocessing.nf`, `conf/modules.config`, `nextflow.config` (add `params.relocate_seed = 42` as the global default).

## Verification

1. `nextflow run . -profile test,docker` succeeds.
2. Outputs under `results/`:
   - N× `<source>.intergenic.bed` (one per source species, no overlap with that source's gene intervals)
   - 2N× `<source>.<ftype>.decoy.gff3`
   - 2N× `<source>.<ftype>.decoy.spliced.fasta`
   - 4N× `<source>.<ftype>[.decoy].spliced.renamed.fasta` (source + decoy variants)
3. Sanity checks:
   - For each source S: `bedtools intersect -a S.intergenic.bed -b S.genes.bed -u | wc -l` → 0
   - For each decoy GFF3: every transcript's `seqid` is a chromosome of the **same source species** (not target, not another source)
   - Decoy spliced FASTA record count matches transcript count in matching decoy GFF3
   - Renamed FASTA record count matches pre-rename record count (no records dropped)
   - Renamed headers match `^>\S+\|(lncRNA|mRNA|decoy_lncRNA|decoy_mRNA)\|<species>$`
   - Re-running with the same `--seed` produces byte-identical decoy GFF3

## Deliberate omissions

- **Step 01** (non-overlapping lncRNA filter) — skipped, consistent with longest-isoform deferral.
- **Step 05** (AGAT statistics) — out of scope; will land in a follow-up alongside MultiQC integration.
- **Multiple targets** — single-target model retained; samplesheet still accepts exactly one `target` row.
- **Strand-aware relocation** — legacy script ignores strand collisions between decoy and existing features; we keep that behavior for parity.
