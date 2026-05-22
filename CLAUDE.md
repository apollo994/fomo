# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Project Context

FOMO is a Nextflow DSL2 pipeline that annotates candidate long non-coding RNAs (lncRNAs) on a target genome assembly by transferring lncRNA annotations from source species. Input: target `.fa` assembly + source `.gff3`/`.fa`. Output: `.gff3` of candidate lncRNAs on the target.

When developing pipelines:

1. **Always use the Seqera MCP tools** to search nf-core modules before writing custom processes.
2. Prefer nf-core module conventions (containers, versioning, meta maps).
3. Use Wave containers for environment management.
4. After writing a pipeline, validate config against the Seqera Platform before launching.
5. Follow DSL2 syntax.

# Nomenclature
- **lncRNA**: long non-coding RNA
- **mRNA**: messenger RNA, genes with CDS annotated
- **target**: assembly/species to be annotated
- **source**: assembly/species with annotation to be transferred to the target

# Pipeline Architecture

The pipeline has five stages (only **preprocessing** is currently being implemented):

1. **Collect source annotation** — gather lncRNA/mRNA GFF3 + FASTA from source species
2. **Preprocessing** — prepare spliced sequences and decoy sequences for alignment
3. **Projection** — minimap2 splice-aware alignment of source sequences onto target
4. **Validation** — splice-junction validation and gffcompare benchmarking
5. **Reporting** — statistics and plots

### Preprocessing detail (current focus)

Each step maps to a legacy script in `legacy_scripts/preprocessing/` which serves as the reference implementation:

| Step | Legacy script | Tool | Purpose |
|------|--------------|------|---------|
| Extract lncRNA/mRNA spliced FASTA | `00_get_feature_annotation.sh` | AGAT | Filter GFF3 by feature type, keep longest isoforms, extract exon sequences |
| Extract intergenic intervals | `02_get_intergenic_intervals.sh` | AGAT + bedtools | Build intergenic BED from target annotation |
| Build decoy sequences | `03_relocate_loci.py` | Python | Relocate source loci into intergenic intervals of target to create decoy FASTA |
| Extract decoy spliced FASTA | `04_get_decoy_sequence_commands.sh` | AGAT | Extract exon sequences from decoy GFF3 |
| Preprocessing statistics | `05_get_gff_statistics_commands.sh` | AGAT | Collect GFF stats at each stage |

The projection stage uses minimap2 (see `legacy_scripts/minimap_transfer/`) and converts BAM → GFF3 with exon structure.

# Key Tools
- **gffread** - GFF3 manipulation. Prefer this over **AGAT** when possible. 
- **AGAT** — GFF3 manipulation (filtering, longest isoform selection, statistics).
- **minimap2** — splice-aware long-read alignment (source spliced FASTA → target)
- **bedtools** — genomic interval arithmetic
- **samtools** — BAM handling
- **gffcompare** — annotation comparison for benchmarking
- **
