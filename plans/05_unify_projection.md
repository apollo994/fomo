# Plan: Unify MAPPING + BAM_TO_GFF under PROJECTION + record alignment %

## Context

Two related changes on top of `plans/03_mapping.md` + `plans/04_bam_to_gff.md`:

1. **Subworkflow restructure.** `BRAINSTORM.md` lists "projection" as a single stage. Today the implementation splits it into two subworkflows:
   - `MAPPING` — `GUNZIP_TARGET → MINIMAP2_INDEX → MINIMAP2_ALIGN`
   - `PROJECTION` — `BAM_TO_GFF`

   Merge into a single `PROJECTION` subworkflow that ingests preprocessed source FASTAs + target FASTA and emits projected GFF3s. Keeps the pipeline shape aligned with the conceptual stages.

2. **Length/coverage attributes on transcript rows.** The current GFF carries `NM/AS/de` (raw aligner tags) but nothing about how much of the source transcript actually projected onto the target. Add `original_length`, `projected_length`, and `pct_align` to every `transcript` row.

These changes ship together because they touch the same files (one process rename + one script change).

## Change 1 — Unify subworkflows

### Decisions

| Decision | Choice |
|----------|--------|
| Keep `MAPPING` subworkflow? | **No.** Merge all four processes under `PROJECTION`. |
| New process namespace | `FOMO:PROJECTION:{GUNZIP_TARGET, MINIMAP2_INDEX, MINIMAP2_ALIGN, BAM_TO_GFF}` |
| Delete `subworkflows/local/mapping.nf`? | **Yes** — replaced entirely, no callers besides `workflows/fomo.nf`. |
| Plan history | Leave `plans/03_mapping.md` + `plans/04_bam_to_gff.md` in place as the original record. This plan supersedes the structural decision they made about subworkflow split. |

### Files touched

- `subworkflows/local/projection.nf` — expand from 1-process wrapper to full chain.
- `subworkflows/local/mapping.nf` — delete.
- `workflows/fomo.nf` — remove `include { MAPPING }`, drop `MAPPING(...)` call, change `PROJECTION(MAPPING.out.bam)` to `PROJECTION(ch_input.target, PREPROCESSING.out.spliced_fasta, PREPROCESSING.out.decoy_spliced_fasta)`.
- `conf/modules.config` — rename four `FOMO:MAPPING:*` selectors to `FOMO:PROJECTION:*` (no behavior change beyond the prefix; same `ext.args`, same `ext.prefix`).

### New `subworkflows/local/projection.nf`

```groovy
include { GUNZIP as GUNZIP_TARGET } from '../../modules/nf-core/gunzip/main'
include { MINIMAP2_INDEX          } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN          } from '../../modules/nf-core/minimap2/align/main'
include { BAM_TO_GFF              } from '../../modules/local/bam_to_gff'

workflow PROJECTION {
    take:
    ch_target            // [ meta(role:'target'), fasta.gz, gff3.gz ] × 1
    ch_spliced_fasta     // [ meta(decoy:false, feature_type), fasta ] × 2N
    ch_decoy_spliced     // [ meta(decoy:true,  feature_type), fasta ] × 2N

    main:

    GUNZIP_TARGET(
        ch_target.map { meta, fa, _gff3 -> tuple(meta, fa) }
    )

    MINIMAP2_INDEX(GUNZIP_TARGET.out.gunzip)

    ch_reference = MINIMAP2_INDEX.out.index
        .map { meta, mmi -> tuple([id: meta.id], mmi) }
        .first()

    ch_reads = ch_spliced_fasta
        .mix(ch_decoy_spliced)
        .combine(ch_reference)
        .map { meta, reads, ref_meta, _mmi ->
            tuple(meta + [target_id: ref_meta.id], reads)
        }

    MINIMAP2_ALIGN(
        ch_reads,
        ch_reference,
        true,    // bam_format
        'bai',   // bam_index_extension
        false,   // cigar_paf_format
        false    // cigar_bam
    )

    BAM_TO_GFF(MINIMAP2_ALIGN.out.bam)

    emit:
    bam   = MINIMAP2_ALIGN.out.bam    // [ meta, *.bam     ]
    index = MINIMAP2_ALIGN.out.index  // [ meta, *.bam.bai ]
    gff3  = BAM_TO_GFF.out.gff3       // [ meta, *.gff3    ]
}
```

### `workflows/fomo.nf` diff (logical)

```groovy
// remove
include { MAPPING } from '../subworkflows/local/mapping'
...
MAPPING(
    ch_input.target,
    PREPROCESSING.out.spliced_fasta,
    PREPROCESSING.out.decoy_spliced_fasta
)
PROJECTION(MAPPING.out.bam)

// becomes
PROJECTION(
    ch_input.target,
    PREPROCESSING.out.spliced_fasta,
    PREPROCESSING.out.decoy_spliced_fasta
)
```

### `conf/modules.config` rename

All four `FOMO:MAPPING:*` selectors → `FOMO:PROJECTION:*`. No other change. The existing `FOMO:PROJECTION:BAM_TO_GFF` block stays as is.

---

## Change 2 — Record original vs projected length on transcript rows

### Decisions

| Decision | Choice |
|----------|--------|
| `original_length` | Sum of CIGAR query-consuming ops (`M, I, S, =, X`) — i.e. the source read length. For spliced FASTA reads this equals the spliced transcript length. Computed from CIGAR (more robust than `length($10)` which can be `*` on bizarre records; `-F 2308` already excludes the usual culprits). |
| `projected_length` | Sum of exon block lengths on the target reference — i.e. the sum of `M/=/X/D` ops that fall inside genomic blocks (already accumulated in the current CIGAR walk; equal to ∑ block_end − block_start + 1). Excludes intronic `N` gaps. |
| `pct_align` | `100 × projected_length / original_length`, formatted to 2 decimal places (`%.2f`). Sentinel `.` if `original_length == 0` (defensive — shouldn't happen). |
| Row scope | `transcript` only. Not on `gene` or `exon` — gene span ≡ transcript span here (one-transcript-per-gene), and exons inherit context via `Parent`. |
| Attribute order on `transcript` row | `ID;Parent;source;gene_class;original_length;projected_length;pct_align;NM;AS;de` — length/coverage attributes before raw aligner tags. |

### `bin/bam_to_gff.sh` changes

Two small additions inside the CIGAR walk:

```awk
query_len     = 0      # sum of M, I, S, =, X
projected_len = 0      # sum across blocks
...
if (op == "M" || op == "=" || op == "X") {
    if (block_len == 0) block_start = ref_pos
    block_len     += len
    ref_pos       += len
    query_len     += len
    projected_len += len
}
else if (op == "D") {
    if (block_len == 0) block_start = ref_pos
    block_len     += len
    ref_pos       += len
    projected_len += len      # D consumes ref but not query
}
else if (op == "N") { ... }              # unchanged
else if (op == "I") { query_len += len } # consumes query only
else if (op == "S") { query_len += len } # soft-clip, consumes query only
else if (op == "H" || op == "P") { }     # unchanged (consume neither)
```

Then on the transcript row:

```awk
pct_align = (query_len > 0 ? sprintf("%.2f", 100.0 * projected_len / query_len) : ".")
tx_attrs  = "ID=" tid ";Parent=gene-" tid \
            ";source=" src ";gene_class=" gclass \
            ";original_length=" query_len \
            ";projected_length=" projected_len \
            ";pct_align=" pct_align \
            ";NM=" nm ";AS=" as ";de=" de
```

Gene + exon attributes unchanged.

### Sample transcript row after the change

Before:
```
HG995163.1  fomo  transcript  21417  49107  41  -  .  ID=ENSGXNT00000003663;Parent=gene-ENSGXNT00000003663;source=Lycaena_alciphron_282377;gene_class=mRNA;NM=146;AS=137;de=0.2202
```

After:
```
HG995163.1  fomo  transcript  21417  49107  41  -  .  ID=ENSGXNT00000003663;Parent=gene-ENSGXNT00000003663;source=Lycaena_alciphron_282377;gene_class=mRNA;original_length=783;projected_length=592;pct_align=75.61;NM=146;AS=137;de=0.2202
```

---

## Verification

1. `nextflow run . -profile test,docker -resume`. Cache should hold for all upstream PREPROCESSING + STATS tasks; `GUNZIP_TARGET` / `MINIMAP2_INDEX` / `MINIMAP2_ALIGN` should also cache (same process scripts, just under a renamed subworkflow). `BAM_TO_GFF` re-runs because the script changed.
2. Channel inspection: `results/bam_to_gff/` still contains 16 `.projected.gff3` files; counts of gene/transcript/exon rows unchanged.
3. Spot check one transcript row contains `original_length`, `projected_length`, `pct_align` and the values look sane:
   - `0 < pct_align ≤ 100` for the vast majority of mRNA primaries
   - `original_length` matches `length(SEQ)` from `samtools view -F 2308 <bam> | awk -v tid=... '$1 ~ tid { print length($10); exit }'` for at least one transcript
4. No process namespaced `FOMO:MAPPING:*` appears in the Nextflow execution report.

## Out of scope

- Filtering by `pct_align` threshold — separate downstream step.
- Adding length/coverage to the gene row — gene span equals transcript span here, would be redundant.
- Changes to alignment QC / samtools stats — still deferred.
