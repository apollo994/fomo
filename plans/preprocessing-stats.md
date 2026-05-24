# Plan: FOMO Preprocessing — Stats & HTML reporting

## Context

Preprocessing now emits filtered GFFs, decoy GFFs, source spliced FASTAs and decoy spliced FASTAs. There's no QC pass over the inputs or outputs, so it's impossible to tell — quickly and visually — whether a run produced sensible annotations or whether something silently degraded (e.g. a source assembly missing lncRNA records, or a decoy pool collapsing to zero on a chromosome). This plan adds a STATS subworkflow that runs `agat_sp_statistics` on every GFF that matters (raw inputs + preprocessing outputs), FastQC on every spliced transcript FASTA (source + decoy), and aggregates everything into a single MultiQC HTML report. Output: one navigable HTML per run that lets the user audit inputs against outputs in one place.

## Decisions (user-confirmed)

| Decision | Choice |
|----------|--------|
| GFF stats scope | All GFFs: raw source (N) + raw target (1) + filtered (2N) + decoy (2N) = **5N+1** AGAT tasks |
| FASTA QC | FastQC on source spliced + decoy spliced = **4N** FASTQC tasks |
| HTML aggregator | MultiQC; FastQC ingested natively, AGAT via `custom_content` |
| AGAT output format | **YAML** (`agat_sp_statistics.pl --yaml`) — structured, easy to `yaml.safe_load()` in the parser; far more robust than the multi-section indented TXT format |
| AGAT module shape | **Local module** (`modules/local/agat_spstatistics.nf`) — the nf-core module captures only `*.txt`, so we fork to a 25-line local wrapper that emits both `.txt` (for human inspection) and `.yaml` (for the parser). Reuses the same AGAT biocontainer that nf-core uses. |
| AGAT → MultiQC glue | Python parser `bin/agat_to_mqc.py` (reads YAML, emits `*_mqc.tsv`) |
| Channel topology for AGAT | Single `AGAT_SPSTATISTICS` process, mixed channel of all 5N+1 GFFs; meta carries `kind: 'raw' \| 'filtered' \| 'decoy'` to drive `ext.prefix` |
| Subworkflow placement | New `subworkflows/local/stats.nf`; invoked from `workflows/fomo.nf` after `PREPROCESSING` |
| Raw GFF decompression | Alias the existing nf-core gunzip module as `GUNZIP_RAW_GFF` — AGAT does not auto-decompress `.gff3.gz` |

## Channel topology

```
ch_input.source (N) ──┐
ch_input.target (1) ──┴─► map: [meta + kind:'raw', gff3.gz]
                           │
                           └─ GUNZIP_RAW_GFF ──► [meta(raw), gff3]   (N+1)

preprocessing.filtered_gff3 ──► map: [meta + kind:'filtered', gff3]  (2N)
preprocessing.decoy_gff3    ──► map: [meta + kind:'decoy',    gff3]  (2N)

   ├──── mix ────► AGAT_SPSTATISTICS ──► [meta, *.stats.txt + *.stats.yaml]   (5N+1)
                                            │
                                            └─ AGAT_TO_MQC(yaml) ──► [meta, *_mqc.tsv]

preprocessing.spliced_fasta       ──┐
preprocessing.decoy_spliced_fasta ──┴─► mix ──► FASTQC ──► [meta, *.zip, *.html]   (4N)

   FASTQC.out.zip + AGAT_TO_MQC.out.mqc ──► collect ──► MULTIQC ──► multiqc_report.html
```

## Files to add

```
subworkflows/local/stats.nf                              # new subworkflow
modules/local/agat_spstatistics.nf                       # local AGAT wrapper emitting txt + yaml
modules/nf-core/fastqc/{main.nf, environment.yml}        # fetched from nf-core/modules master
modules/nf-core/multiqc/{main.nf, environment.yml}       # fetched
modules/local/agat_to_mqc.nf                             # thin wrapper around bin/agat_to_mqc.py
bin/agat_to_mqc.py                                       # AGAT YAML → MultiQC _mqc.tsv parser
assets/multiqc_config.yml                                # report title, custom_content section names, sample-name cleaning rules
```

## Files to edit

```
workflows/fomo.nf                # wire STATS after PREPROCESSING, pass ch_input.source + ch_input.target
conf/modules.config              # withName entries for AGAT_SPSTATISTICS, FASTQC, MULTIQC, AGAT_TO_MQC, GUNZIP_RAW_GFF
nextflow.config                  # declare params.multiqc_config default = "${projectDir}/assets/multiqc_config.yml"
```

## Local AGAT module (`modules/local/agat_spstatistics.nf`)

```groovy
process AGAT_SPSTATISTICS {
    tag "${meta.id}.${meta.kind}"
    label 'process_single'

    container 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.stats.txt"),  emit: stats_txt
    tuple val(meta), path("*.stats.yaml"), emit: stats_yaml
    tuple val("${task.process}"), val('agat'),
          eval('agat_sp_statistics.pl --help | grep -i version | head -n1 | sed "s/[^0-9.]//g"'),
          topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    agat_sp_statistics.pl \\
        --gff ${gff} \\
        --output ${prefix}.stats.txt \\
        --yaml ${prefix}.stats.yaml \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    touch ${prefix}.stats.txt ${prefix}.stats.yaml
    """
}
```

## Subworkflow signature (`subworkflows/local/stats.nf`)

```groovy
include { GUNZIP as GUNZIP_RAW_GFF } from '../../modules/nf-core/gunzip/main'
include { AGAT_SPSTATISTICS        } from '../../modules/local/agat_spstatistics'
include { AGAT_TO_MQC              } from '../../modules/local/agat_to_mqc'
include { FASTQC                   } from '../../modules/nf-core/fastqc/main'
include { MULTIQC                  } from '../../modules/nf-core/multiqc/main'

workflow STATS {
    take:
    ch_raw_source          // [meta, fasta.gz, gff3.gz] × N
    ch_raw_target          // [meta, fasta.gz, gff3.gz] × 1
    ch_filtered_gff3       // [meta(decoy:false), gff3] × 2N
    ch_decoy_gff3          // [meta(decoy:true),  gff3] × 2N
    ch_spliced_fasta       // [meta(decoy:false), fasta] × 2N
    ch_decoy_spliced_fasta // [meta(decoy:true),  fasta] × 2N

    main:
    // Decompress raw GFFs (source + target) — AGAT can't read .gz
    ch_raw_gff_gz = ch_raw_source.mix(ch_raw_target)
        .map { meta, fa, gff -> tuple(meta + [kind: 'raw'], gff) }
    GUNZIP_RAW_GFF(ch_raw_gff_gz)

    // Mix raw + filtered + decoy → AGAT
    ch_agat_in = GUNZIP_RAW_GFF.out.gunzip
        .mix(ch_filtered_gff3.map { m, g -> tuple(m + [kind: 'filtered'], g) })
        .mix(ch_decoy_gff3   .map { m, g -> tuple(m + [kind: 'decoy'],    g) })
    AGAT_SPSTATISTICS(ch_agat_in)

    // Parse YAML output into MultiQC custom-content TSV
    AGAT_TO_MQC(AGAT_SPSTATISTICS.out.stats_yaml)

    // FastQC on all spliced FASTAs
    FASTQC(ch_spliced_fasta.mix(ch_decoy_spliced_fasta))

    // MultiQC — collect everything
    ch_mqc = FASTQC.out.zip.map { _m, z -> z }
        .mix(AGAT_TO_MQC.out.mqc.map { _m, t -> t })
        .collect()
    MULTIQC(
        ch_mqc,
        file(params.multiqc_config, checkIfExists: true),
        [], [], []
    )

    emit:
    agat_stats_txt  = AGAT_SPSTATISTICS.out.stats_txt
    agat_stats_yaml = AGAT_SPSTATISTICS.out.stats_yaml
    fastqc_html     = FASTQC.out.html
    multiqc_html    = MULTIQC.out.report
}
```

## `conf/modules.config` additions

```groovy
withName: 'FOMO:STATS:GUNZIP_RAW_GFF' {
    ext.prefix = { "${meta.id}.raw" }
}
withName: 'FOMO:STATS:AGAT_SPSTATISTICS' {
    ext.prefix = { [
        meta.id,
        meta.kind == 'raw'      ? 'raw'
      : meta.kind == 'filtered' ? meta.feature_type
      :                           "${meta.feature_type}.decoy"
    ].join('.') }
    // Raw whole-genome GFFs are much larger than filtered/decoy
    memory = { meta.kind == 'raw' ? 8.GB : 2.GB }
    time   = { meta.kind == 'raw' ? 1.h  : 15.min }
}
withName: 'FOMO:STATS:AGAT_TO_MQC' {
    ext.prefix = { [
        meta.id,
        meta.kind == 'raw'      ? 'raw'
      : meta.kind == 'filtered' ? meta.feature_type
      :                           "${meta.feature_type}.decoy"
    ].join('.') }
}
withName: 'FOMO:STATS:FASTQC' {
    ext.prefix = { "${meta.id}.${meta.feature_type}${meta.decoy ? '.decoy' : ''}.spliced" }
}
```

## `bin/agat_to_mqc.py` — YAML parser

`agat_sp_statistics.pl --yaml` produces a structured YAML doc. The exact schema (verified against the legacy outputs in the project) is roughly:

```yaml
mrna:
  number of mrnas: 1234
  number of mrnas with utr: 567
  mean mrna length (bp): 1200.5
  total mrna length (bp): 1481234
exon:
  number of exons: 5678
  ...
gene:
  number of genes: 1000
  ...
```

Parser logic:

1. `yaml.safe_load(open(args.input))` → nested dict
2. Flatten into a single sample row: keys = `<section>__<metric>` (e.g. `mrna__number_of_mrnas`); whitespace replaced with `_`.
3. Cherry-pick a curated subset of metrics for the MultiQC table (full set is overwhelming): `gene__number_of_genes`, `mrna__number_of_mrnas`, `exon__number_of_exons`, `mrna__mean_mrna_length_bp`, `mrna__total_mrna_length_bp`, `cds__number_of_cdss`. Remaining metrics retained in a secondary `_all_mqc.tsv` for users who want depth.
4. Write a `*_mqc.tsv` file with the MultiQC custom-content header block:
   ```
   # id: agat_summary
   # section_name: 'AGAT annotation statistics'
   # description: 'Per-annotation summary from agat_sp_statistics (--yaml)'
   # plot_type: 'table'
   # pconfig:
   #     id: 'agat_summary_table'
   #     namespace: 'AGAT'
   Sample\tn_gene\tn_mrna\tn_exon\tn_cds\tmean_mrna_length\ttotal_mrna_length
   <name>\t...
   ```

CLI:
```
agat_to_mqc.py --input <stats.yaml> --name <prefix> --output <prefix>_mqc.tsv
```

If the YAML schema differs across AGAT versions, the parser logs a warning and emits whatever flat metrics it found — the test profile will catch regressions via the verification step.

The wrapper module `modules/local/agat_to_mqc.nf` is a ~25-line process (python:3.11-slim container, label `process_single`, ships `pyyaml` via a one-line `pip install pyyaml` in the script block — or, cleaner, use a community Wave container that already has pyyaml. Default to the inline `pip install` for first pass since it's a 3 MB dep).

## `assets/multiqc_config.yml` — presentation

```yaml
report_title: "FOMO preprocessing report"
report_comment: "Per-source/target GFF stats (AGAT) + transcript FASTA QC (FastQC)"
extra_fn_clean_exts:
  - ".spliced.renamed"
  - ".spliced"
  - ".decoy"
  - ".raw"
  - "_mqc"
custom_data:
  agat_summary:
    file_format: "tsv"
    section_name: "AGAT annotation statistics"
    plot_type: "table"
sp:
  agat_summary:
    fn: "*agat_summary_mqc.tsv"
fastqc_config:
  fastqc_theoretical_gc: ""  # disable GC reference panel — meaningless for transcript FASTAs
```

## Workflow wiring (`workflows/fomo.nf`)

Add an assertion that exactly one target row exists, then invoke STATS:

```groovy
include { STATS } from '../subworkflows/local/stats'

// after PREPROCESSING(...)
STATS(
    ch_input.source,
    ch_input.target,
    PREPROCESSING.out.filtered_gff3,
    PREPROCESSING.out.decoy_gff3,
    PREPROCESSING.out.spliced_fasta,
    PREPROCESSING.out.decoy_spliced_fasta
)
```

## Reused utilities

- `modules/nf-core/gunzip/main.nf` — re-aliased as `GUNZIP_RAW_GFF`, no new install needed
- The `meta + [kind: 'X']` pattern mirrors `meta + [feature_type: ftype]` and `meta + [decoy: true]` already used in `subworkflows/local/preprocessing.nf`
- The `.mix(...)` + meta-discriminated single-process pattern follows `RENAME_FASTA_HEADERS` (preprocessing subworkflow)
- AGAT biocontainer (`quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1`) — matches the nf-core/agat/spstatistics version

## Verification

1. `nextflow run main.nf -profile test,docker` succeeds end-to-end.
2. `results/agat_spstatistics/` contains **2 × (5N+1) = 42** files for N=4 test sources — 21 `.stats.txt` + 21 `.stats.yaml`.
3. Spot-check one YAML: `python3 -c 'import yaml; print(list(yaml.safe_load(open("...stats.yaml")).keys()))'` returns a list including `gene`, `mrna`, `exon`.
4. `results/agat_to_mqc/` contains **21** `*_mqc.tsv` files; each has the MultiQC header block + one data row.
5. `results/fastqc/` contains **4N = 16** .html + .zip pairs.
6. `results/multiqc/multiqc_report.html` exists. Opening it shows:
   - A FastQC section with 16 samples
   - An "AGAT annotation statistics" custom table with 21 rows
   - Sample names cleaned of `.spliced.renamed` / `.spliced` / `.decoy` / `.raw` suffixes
7. Sanity check the AGAT table: filtered `n_mrna` per source should be ≤ raw `n_mrna`; decoy `n_mrna` should equal filtered `n_mrna` (decoys are 1:1 with source transcripts that successfully relocated).
8. Re-run with `-resume` — STATS tasks cache cleanly.

## Deliberate omissions

- **SEQKIT_STATS** — FastQC on FASTA emits many N/A panels; SEQKIT_STATS would give richer length/N50/GC for transcriptome FASTAs. Not in user's explicit ask; can be added later as a second MultiQC custom section.
- **CUSTOM_DUMPSOFTWAREVERSIONS / unified versions.yml** — preprocessing emits topic-based `versions` tuples but doesn't aggregate them into a single `versions.yml` for MultiQC's software-versions panel. Defer to a follow-up.
- **Genome FASTA passed to AGAT** — legacy invoked `agat_sp_statistics.pl -g <fasta>` for chromosome-level stats. Skipped here for simplicity; can be added by extending the local AGAT module's input tuple to optionally include FASTA.
- **MultiQC versions / software panel** — depends on version aggregation above.
