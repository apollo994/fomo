# Plan: FOMO BAM → GFF3 conversion

## Context

`MAPPING` (plan `03_mapping.md`) emits one sorted+indexed BAM per (source × feature_type × decoy?). The next step from `BRAINSTORM.md` is **projection** — converting each BAM into a GFF3 describing the projected transcripts on the target genome. Reference implementation: `legacy_scripts/minimap_transfer/01_convert_bam_to_gff.sh` (driven by `01_convert_bam_to_gff_commands.sh`, which is just a fanout that runs the script on every per-(target, source, class) BAM).

This plan adapts the legacy approach so the resulting GFF3 carries `source` and `gene_class` as first-class attributes — derived by parsing the BAM read names — rather than burying them inside the ID string.

## Read name format

Read names in every BAM come from `RENAME_FASTA_HEADERS` (`modules/local/rename_fasta_headers.nf`) and look like:

```
ENSGXNT00000003939|mRNA|Lycaena_alciphron_282377
ENSFOJT00000002660.reloc|decoy_lncRNA|Lycaena_virgaureae_282395
```

That is: `<tid>|<gene_class>|<source_species>`, where:
- `gene_class ∈ {mRNA, lncRNA, decoy_mRNA, decoy_lncRNA}`
- `source_species` matches `meta.id` upstream
- decoy `tid`s carry a `.reloc` suffix (preserved — keeps the link to `RELOCATE_LOCI`)

## Decisions

| Decision | Choice |
|----------|--------|
| Transcript/gene ID | Plain `<tid>` (with `.reloc` preserved on decoys). Source-and-class are attributes, not part of the ID. |
| Decoy suffix | Keep `.reloc` in the tid — self-documenting link back to `RELOCATE_LOCI`. |
| Exon attributes | Just `ID` + `Parent` (match legacy). Context inherited via `Parent`. |
| Tool versions | Match legacy: `samtools view -F 2308`, `bedtools bamtobed -bed12`, then awk. |
| Filtering | `-F 2308` (drop unmapped, secondary, supplementary) — same as legacy. |
| Output | One `.gff3` per input BAM, with `##gff-version 3` header. |

## Output GFF3 attribute layout

| Row | 9th-column attributes |
|-----|-----------------------|
| `gene` | `ID=gene-<tid>;source=<src>;gene_class=<gclass>` |
| `transcript` | `ID=<tid>;Parent=gene-<tid>;source=<src>;gene_class=<gclass>;NM=<nm>;AS=<as>;de=<de>` |
| `exon` | `ID=<tid>.exonN;Parent=<tid>` |

Other columns unchanged from legacy:
- col 2 (source) = `bedtools` (could later be `fomo` — not changing here)
- col 5 (score) = MAPQ (from BED12 score column)
- col 6 (strand) = from BED12
- coords on `gene` + `transcript` = full alignment span; on `exon` = per-block from BED12 `blockSizes`/`blockStarts`. Zero-size blocks skipped (CIGAR `I` between `N`s).

## 1. New local module: `BAM_TO_GFF`

`modules/local/bam_to_gff.nf` — wraps the legacy two-pass awk pipeline with the new read-name parsing. Conda env supplies `samtools` + `bedtools`. Use a mulled biocontainer (same one used by other modules) or the standard `samtools` + `bedtools` biocontainer combo.

```groovy
process BAM_TO_GFF {
    tag "${meta.target_id}.from_${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/samtools_bedtools:<tag>"   // pick at impl time

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.gff3"), emit: gff3

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam_to_gff.sh ${bam} > ${prefix}.gff3
    """
}
```

The actual logic goes into `bin/bam_to_gff.sh` (Nextflow auto-adds `bin/` to `\$PATH`) — keeps the process body small and the script testable in isolation against any BAM.

## 2. `bin/bam_to_gff.sh`

Two-pass awk like the legacy. Differences vs `legacy_scripts/minimap_transfer/01_convert_bam_to_gff.sh`:

1. Split `$4` (BED12 name / SAM QNAME) on `|` once at the top of the awk block:
   ```awk
   n_parts = split($4, parts, "|")
   tid    = parts[1]
   gclass = (n_parts >= 2 ? parts[2] : ".")
   src    = (n_parts >= 3 ? parts[3] : ".")
   ```
2. Emit IDs as `tid` (not full read name).
3. Append `;source=<src>;gene_class=<gclass>` to gene + transcript rows.
4. Exon rows keep just `ID=<tid>.exonN;Parent=<tid>`.
5. Reject reads whose `tid` (after split) contains `;` or `=` — would corrupt attribute parsing. Legacy already rejects `;` in full name; we tighten to also reject `=`.
6. Second awk pass (NM/AS/de tag injection) keys on `tid` instead of full QNAME — keys are derived from the SAM `QNAME` field by the same split.

## 3. Subworkflow wiring

Extend `subworkflows/local/mapping.nf` (rename to keep the file name, or add a new subworkflow `projection.nf`?). Recommendation: **keep `MAPPING` focused on read-aligner output (BAM) and add a new subworkflow `subworkflows/local/projection.nf`** that takes `MAPPING.out.bam` and runs `BAM_TO_GFF`. Keeps responsibilities separable and matches the naming used in `BRAINSTORM.md` (projection).

```groovy
// subworkflows/local/projection.nf
include { BAM_TO_GFF } from '../../modules/local/bam_to_gff'

workflow PROJECTION {
    take:
    ch_bam   // [ meta(target_id, id, feature_type, decoy), bam ] × 16

    main:
    BAM_TO_GFF(ch_bam)

    emit:
    gff3 = BAM_TO_GFF.out.gff3   // [ meta, *.gff3 ]
}
```

Top-level wiring in `workflows/fomo.nf`:

```groovy
include { PROJECTION } from '../subworkflows/local/projection'
...
MAPPING(...)
PROJECTION(MAPPING.out.bam)
```

## 4. `conf/modules.config` additions

```groovy
withName: 'FOMO:PROJECTION:BAM_TO_GFF' {
    ext.prefix = { "${meta.target_id}.from_${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.projected" }
}
```

## 5. Expected outputs

For the current test samplesheet: **16 GFF3 files** under `results/bam_to_gff/`, one per BAM:

```
Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.mRNA.projected.gff3
Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.mRNA.decoy.projected.gff3
Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.lnc_RNA.projected.gff3
Lycaena_phlaeas_282391.from_Lycaena_alciphron_282377.lnc_RNA.decoy.projected.gff3
... (× 4 sources)
```

Sample transcript row (e.g. mRNA hit):
```
SUPER_1  bedtools  transcript  12345  67890  60  +  .  ID=ENSGXNT00000003939;Parent=gene-ENSGXNT00000003939;source=Lycaena_alciphron_282377;gene_class=mRNA;NM=42;AS=1234;de=0.0123
```

## 6. Verification

1. `nextflow run . -profile test,docker -resume` — confirm `PROJECTION` runs after `MAPPING` and produces 16 GFF3 files.
2. Spot-check one GFF3:
   ```
   awk '$3=="transcript"' results/bam_to_gff/<file>.gff3 | head
   ```
   Confirm `source=...`, `gene_class=...`, and `NM/AS/de` are all present.
3. `grep -c '^[^#]' results/bam_to_gff/<file>.gff3` matches `3 × transcript_count + n_exons` (gene + transcript + exons per read).
4. Diff against legacy on the same input — modulo the ID/attribute rewrite, gene/transcript/exon counts and coordinates should match exactly.

## 7. Deliberate omissions

- **Validation / splice-junction filtering** — next pipeline stage, separate plan.
- **gffcompare benchmarking** — reporting stage.
- **Merging per-source GFF3s into a single annotation** — out of scope; downstream concern.
