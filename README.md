# FOMO

A Nextflow DSL2 pipeline that annotates candidate long non-coding RNAs (lncRNAs) on a
target genome assembly by transferring lncRNA annotations from source species.

- **Input:** target `.fa` assembly + source `.gff3`/`.fa`
- **Output:** `.gff3` of candidate lncRNAs on the target

Each species in the samplesheet can act as a `source` (donates annotation), a `target`
(receives projections), or `both` — so a single run can cover an all-vs-all comparison
across many species.

Pipeline stages:

1. **Preprocessing** — extract spliced lncRNA (and optionally mRNA) sequences from each
   source, plus optional decoy sequences for a false-positive baseline.
2. **Projection** — splice-aware alignment (minimap2) of source sequences onto each
   target genome, converted back to GFF3.
3. **Benchmarking** — compare projected models against each target's existing annotation
   (gffcompare), when available, to score accuracy per source.
4. **Consensus** — rank sources by accuracy and build a top-3 consensus annotation per
   target.
5. **Reporting** — one MultiQC report per target, plus one run-level summary report.

## Usage

```bash
nextflow run . -profile test
```

See `nextflow run . --help` for the full parameter list, and `assets/schema_input.json`
for the samplesheet format (`species,role,fasta,gff3`).

## Documentation

See `CLAUDE.md` for pipeline architecture, subworkflow details, and conventions.
