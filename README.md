# fomo

A Nextflow DSL2 pipeline that assesses genome annotation completeness by running [BUSCO v6](https://busco.ezlab.org/) in transcriptome mode on mRNA FASTA files. Targets Lepidoptera species using the `lepidoptera_odb12` lineage by default.

## Usage

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

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `assets/input/*.mrna.fa` | Glob pattern for input mRNA FASTA files |
| `--outdir` | `results/` | Output directory |
| `--lineage` | `lepidoptera_odb12` | BUSCO lineage database |
| `--busco_mode` | `transcriptome` | BUSCO mode (`genome` / `transcriptome` / `proteins`) |
| `--busco_lineages_path` | `null` | Path to a pre-downloaded lineage database (skips download step) |

## Pipeline structure

```
assets/input/*.mrna.fa
        │
        ├─► BUSCO_DOWNLOAD (lineage DB) ─────────┐
        │                                         │
        └─► BUSCO_BUSCO ◄────────────────────────┘
                │
                └─► results/busco/short_summary.*.json
```

**Modules:**
- `modules/local/busco_download/` — fetches the BUSCO lineage dataset (1 CPU, 16 GB, 1h)
- `modules/nf-core/busco/busco/` — runs BUSCO on each input FASTA in parallel (4 CPUs, 4 GB, 20 min per attempt)

## Outputs

BUSCO JSON short summaries are saved to `results/busco/`. Full tables, sequences, and logs are produced in the working directory but not published by default.

## Module management

nf-core modules are tracked in `modules.json`. To update:

```bash
nf-core modules update busco/busco
nf-core modules install <module-name>
```
