# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A Nextflow DSL 2 pipeline that runs BUSCO v6.0.0 in transcriptome mode on mRNA FASTA files to assess genome annotation completeness. Targets Lepidoptera species using the `lepidoptera_odb12` lineage by default.

## Running the Pipeline

```bash
# Default run (uses input pattern from nextflow.config)
nextflow run main.nf

# Custom parameters
nextflow run main.nf --input 'path/to/*.fa' --outdir results --lineage insects_odb10

# Execution profiles
nextflow run main.nf -profile crg     # CRG HPC cluster (Slurm + Apptainer)
nextflow run main.nf -profile docker  # Docker containers
nextflow run main.nf -profile conda   # Conda environment
```

## Key Parameters (nextflow.config)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `*_mRNA_*.mrna.fa` | Glob pattern for input FASTA files |
| `--outdir` | (set in config) | Results output directory |
| `--lineage` | `lepidoptera_odb12` | BUSCO lineage database |
| `--busco_mode` | `transcriptome` | BUSCO mode (genome/transcriptome/proteins) |

## Architecture

**Workflow (`main.nf`):**
1. Reads input FASTA files matching the configured glob pattern
2. Creates `[meta, fasta]` tuples with file basename as the sample ID
3. `BUSCO_DOWNLOAD` — downloads the specified lineage database once
4. `BUSCO_BUSCO` — runs BUSCO on each FASTA in parallel using the downloaded lineage

**Modules:**
- `modules/local/busco_download/` — local module to fetch BUSCO lineage datasets (1 CPU, 16 GB, 1h)
- `modules/nf-core/busco/busco/` — nf-core BUSCO process (4 CPUs, 4 GB, 20 min/attempt); handles Augustus path configuration for containerized environments

**Key data flow:**
```
assets/input/*.mrna.fa → BUSCO_DOWNLOAD (lineage DB) ─┐
                       └→ BUSCO_BUSCO ←────────────────┘
                              ↓
                   batch_summary, short_summaries (TXT/JSON),
                   full_table, missing_buscos, sequences, logs
```

## Module Management

This pipeline uses nf-core modules tracked in `modules.json`. To update or add modules:
```bash
nf-core modules update busco/busco
nf-core modules install <module-name>
```

## CRG Profile Notes

The `crg` profile uses Slurm for job scheduling and Apptainer (Singularity) for containers. Container images are pulled from the nf-core module's `environment.yml` Conda spec or Docker Hub and cached locally.
