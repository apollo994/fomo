# FOMO

FOMO is a Nextflow pipeline that annotates candidate long non-coding RNAs (lncRNAs) on genome
assemblies by transferring the lncRNA annotation of related species. Each source's spliced
transcripts are aligned to every target genome, the projected models are benchmarked against
the target's own annotation where one exists, and the intron chains that several species agree
on become a curated annotation, merged into the target's reference.

## Requirements

- Nextflow ≥ 25.04
- Singularity / Apptainer or Docker (every tool runs in a pinned container)

## Quick start

```bash
nextflow run . -profile test,singularity     # or: -profile test,docker
```

The test profile runs three species, each cut to one chromosome, in a few minutes.

## Samplesheet

```csv
species,role,fasta,gff3
Lycaena_alciphron_282377,both,alciphron.fna.gz,alciphron.gff3.gz
Lycaena_hippothoe_580924,target,hippothoe.fna.gz,hippothoe.gff3.gz
Lycaena_thersamon_265360,source,thersamon.fna.gz,thersamon.gff3.gz
```

One row per species; `species` is unique and matches `^[A-Za-z0-9_-]+$`. FASTA and GFF3 may be
plain or gzipped.

| role | donates its annotation | receives projections | gff3 |
|:---|:---|:---|:---|
| `source` | yes | no | required |
| `target` | no | yes | optional: the curation reference |
| `both` | yes, to every other target | yes, from every other source | required |

A species is never projected onto itself. Only `both` species are benchmarked: the reference is
the species' own filtered annotation, which only a donor has. A `target`'s gff3 is used for
curation only — curated models overlapping its exons are dropped, each curated gene is tagged
with its nearest reference gene, and the curated genes are merged into it; a `target` without a
gff3 is curated without a reference. One run with `both` rows is an all-vs-all comparison.

## Parameters

`nextflow run . --help` lists every parameter. The main ones:

| parameter | default | |
|:---|:---|:---|
| `--input` | | the samplesheet |
| `--outdir` | `results` | |
| `--include_mrna` | `false` | also transfer mRNA, a positive control |
| `--include_decoy` | `false` | also project relocated decoys, a false-positive baseline |
| `--longest_isoform` | `true` | one transcript per gene, the longest spliced isoform |
| `--curate` / `--curate_min_species` | `true` / `2` | the curated annotation, and how many species must share a chain |
| `--td2_psauron_min` | `0.5` | PSAURON score at which a complete ORF makes a lncRNA coding |

Usage on a cluster, resources and reruns: [docs/usage.md](docs/usage.md).

## Outputs

Per target: the alignments, the projected models (pooled and per source), the benchmarking,
the curated annotation and a MultiQC report; per source: the filtered annotation that was
projected; per run: a summary report and its tables. The full tree: [docs/output.md](docs/output.md).

## Documentation

- [docs/usage.md](docs/usage.md) — running FOMO
- [docs/output.md](docs/output.md) — what each output file holds
- [docs/architecture.md](docs/architecture.md) — how the pipeline works
- [docs/data_records.md](docs/data_records.md) — the statistics records behind the reports
- [CONTRIBUTING.md](CONTRIBUTING.md) — conventions, tests and gates
- [CHANGELOG.md](CHANGELOG.md)

## Citation and licence

There is no publication yet; cite this repository. MIT licence — see [LICENSE](LICENSE).
