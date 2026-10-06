# FOMO

A Nextflow DSL2 pipeline that annotates candidate long non-coding RNAs (lncRNAs) on a
target genome assembly by transferring lncRNA annotations from source species.

- **Input:** target `.fa` assembly + source `.gff3`/`.fa`
- **Output:** `.gff3` of candidate lncRNAs on the target

Each species in the samplesheet can act as a `source` (donates annotation), a `target`
(receives projections), or `both` — so a single run can cover an all-vs-all comparison
across many species. A `gff3` is required for `source` and `both` and must be left empty
for `target`: to benchmark a species against its own annotation, make it `both`.

Pipeline stages:

1. **Preprocessing** — keep one spliced transcript per gene (the longest isoform;
   `--longest_isoform false` keeps them all), drop lncRNA with coding potential (TD2), and
   extract their sequences from each source, plus optional decoy sequences for a
   false-positive baseline. The result is one filtered annotation per source.
2. **Projection** — splice-aware alignment (minimap2) of source sequences onto each
   target genome, converted back to GFF3.
3. **Benchmarking** — compare projected models against each `both` target's own filtered
   annotation (gffcompare) to score accuracy per source.
4. **Consensus** — rank sources by accuracy and build a top-3 consensus annotation per
   target.
5. **Curation** — keep only intron chains shared exactly by ≥ 2 source species
   (`--curate_min_species`), drop models overlapping any exon of the target's reference,
   report one representative per gene (tagged with its location relative to the
   reference) and merge the result into the reference annotation (`--curate false` skips it).
6. **Reporting** — one MultiQC report per target, plus one run-level summary report.

## Usage

```bash
nextflow run . -profile test
```

See `nextflow run . --help` for the full parameter list, and `assets/schema_input.json`
for the samplesheet format (`species,role,fasta,gff3`).

## Outputs

```
targets/<target>/   alignment/, annotation/, gffcompare/, select_top_sources/, multiqc/,
                    curated/<target>.curated.<gtype>.gff3.gz          — curated genes
                            <target>.curated.<gtype>.merged.gff3.gz   — reference + curated
                            <target>.curated.<gtype>.curation.{tsv,json} — curation report
sources/<source>/   annotation/<source>.<ft>.filtered.gff3 — what was projected from
                    this source (and, for a `both` species, the benchmark reference)
summary/            multiqc/ (run-level report), tables/
```

## Documentation

See `CLAUDE.md` for pipeline architecture, subworkflow details, and conventions.

## License

MIT — see [LICENSE](LICENSE).
