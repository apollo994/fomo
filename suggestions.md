
The content of this file should serve as inpiration for the pipeline design.



# CLAUDE.md — nf-core/crossmap

## Project Overview

nf-core/crossmap is a Nextflow DSL2 pipeline for mapping genomic features across species. Built from the nf-core/tools template v3.2.0.

- **Nextflow:** >= 24.04.2
- **Plugin:** nf-schema v2.3.0 (parameter validation and samplesheet parsing)
- **License:** MIT

## Quick Commands

```bash
# Run pipeline with test data
nextflow run main.nf -profile docker,test --outdir results

# Run nf-test suite
nf-test test --profile debug,test,docker --verbose

# Lint pipeline against nf-core standards
nf-core pipelines lint

# Install/update an nf-core module
nf-core modules install <tool_name>
nf-core modules update <tool_name>

# Rebuild parameter schema after adding/changing params
nf-core pipelines schema build
```

**Available profiles:** `test`, `test_full`, `docker`, `singularity`, `conda`, `podman`, `shifter`, `charliecloud`, `apptainer`, `wave`

## Architecture

### Entry Point

`main.nf` orchestrates three phases:

1. **PIPELINE_INITIALISATION** — validates parameters, parses samplesheet, creates channels
2. **NFCORE_CROSSMAP** — runs the main analysis workflow (from `workflows/crossmap.nf`)
3. **PIPELINE_COMPLETION** — sends email/Slack notifications, prints summary

### Directory Layout

```
main.nf                     # Entry point
nextflow.config             # Master config (params, profiles, plugins)
nextflow_schema.json        # JSON Schema for parameter validation
modules.json                # Tracks installed nf-core modules/subworkflows with git SHAs
.nf-core.yml                # nf-core template metadata (org, name, version)

workflows/
  crossmap.nf               # Primary workflow logic

modules/
  nf-core/                  # Reusable process modules from nf-core registry
    minimap2/align/         # Spliced alignment (minimap2 -ax splice)
    multiqc/                # Report aggregation
  local/                    # Pipeline-specific process modules
    filter_features.nf      # Extract feature types from GFF3 (awk + AGAT)
    agat_longest_isoform.nf # Keep longest isoform per gene
    agat_extract_sequences.nf # GFF3 + FASTA → exon sequences
    get_non_overlapping.nf  # lncRNA not overlapping mRNA (bedtools + AGAT)
    get_intergenic_intervals.nf # Identify intergenic regions (AGAT)
    relocate_loci.nf        # Decoy generation (STUB — needs user logic)
    bam_to_gff.nf           # BAM → GFF3 conversion (STUB — needs user logic)
    annocli_alias.nf        # Alias normalization (STUB — needs annocli)
    annocli_download.nf     # Data download (STUB — needs annocli)

subworkflows/
  nf-core/                  # Utility subworkflows from nf-core registry
    utils_nextflow_pipeline/    # Version display, parameter dumping
    utils_nfcore_pipeline/      # Config validation, completion handlers
    utils_nfschema_plugin/      # Parameter validation against schema
  local/
    utils_nfcore_crossmap_pipeline/  # Pipeline-specific init/completion logic
    prepare_genome/          # Steps 1-2: download + alias normalization
    extract_features/        # Steps 3-7: filter → longest → non-overlap → decoy → extract

conf/
  base.config               # Default resource allocation (CPU, memory, time per label)
  modules.config             # Per-module ext.args and publishDir settings
  igenomes.config            # Reference genome paths (iGenomes)
  test.config                # Minimal test profile (small data)
  test_full.config           # Full-size test profile

assets/                     # Email/Slack templates, MultiQC config, samplesheet schema, logos
docs/                       # usage.md, output.md
```

## Key Conventions

### Naming

- **Workflows / Subworkflows / Processes:** UPPERCASE with underscores — `NFCORE_CROSSMAP`, `FASTQC`
- **Functions:** camelCase — `validateInputParameters()`, `getGenomeAttribute()`
- **Channels:** `ch_` prefix — `ch_versions`, `ch_samplesheet`, `ch_multiqc_files`
- **Channel naming pattern:** `ch_output_from_<process>` or `ch_<prev>_for_<next>`
- **Files and directories:** snake_case

### Module Structure

Each module in `modules/nf-core/<name>/` contains:

| File | Purpose |
|------|---------|
| `main.nf` | Process definition (container specs, input/output, script, stub) |
| `meta.yml` | Tool metadata, authors, I/O documentation |
| `environment.yml` | Conda environment |
| `tests/` | nf-test specifications |

### Configuration

- **Process resources** — assigned via labels in `conf/base.config`: `process_single` (1 CPU, 6 GB), `process_low` (2 CPU, 12 GB), `process_medium` (6 CPU, 36 GB), `process_high` (12 CPU, 72 GB), `process_high_memory` (200 GB), `process_long` (20h)
- **Module arguments** — set via `ext.args` / `ext.args2` / `ext.args3` in `conf/modules.config`
- **Output publish paths** — set via `publishDir` in `conf/modules.config`
- **Version tracking** — every process must emit a `versions.yml` file

### Configuration Precedence (lowest to highest)

1. `nextflow.config` defaults
2. `conf/base.config` (resource labels)
3. Profile-specific configs (`conf/test.config`, etc.)
4. `conf/igenomes.config`
5. `conf/modules.config` (loaded last in nextflow.config)
6. User custom config (`-c /path/to/config`)
7. Command-line parameters (`--param value`)

## Development Workflow

- Contribute to `dev` branch; patches go to `main`
- Run `nf-core pipelines lint` before submitting PRs
- Run `nf-test test --profile debug,test,docker --verbose` for testing
- Rebuild `nextflow_schema.json` via `nf-core pipelines schema build` when adding parameters
- Update `docs/output.md` when adding new output files
- Update `CITATIONS.md` when adding new tools
- Update `CHANGELOG.md` with changes
- Pre-commit hooks enforce Prettier formatting and EditorConfig compliance

## CI/CD

- **ci.yml** — Matrix testing across Nextflow versions (24.04.2, latest) and container profiles (docker, singularity, conda)
- **linting.yml** — Prettier formatting + nf-core lint checks
- **awstest.yml / awsfulltest.yml** — AWS batch testing via Seqera Platform

## Pipeline-Specific Design

### Samplesheet Format

```csv
sample,taxid,role,assembly,annotation
Pieris_rapae,64459,both,/path/to/genome.fna,/path/to/annotation.gff3
Pieris_brassicae,7116,source,/path/to/genome.fna,/path/to/annotation.gff3
Vanessa_cardui,110448,target,/path/to/genome.fna,
```

- **role:** `source` (provides gene models), `target` (receives mapped models), `both`
- **annotation:** required for source/both roles
- Channel creation splits samplesheet into `ch_sources` and `ch_targets`

### Pipeline Steps

1. ANNOCLI_DOWNLOAD — fetch genome + annotation by taxid (optional)
2. ANNOCLI_ALIAS — normalize sequence identifiers between GFF3 and FASTA
3. FILTER_FEATURES — extract feature types (lnc_RNA, mRNA) from source annotations
4. AGAT_LONGEST_ISOFORM — keep longest isoform per gene
5. GET_NON_OVERLAPPING — lncRNAs not overlapping mRNA genes (bedtools intersect)
6. GENERATE_DECOYS — relocate lncRNAs to intergenic regions as negative controls
7. EXTRACT_SEQUENCES — GFF3 + genome FASTA → exon FASTA (AGAT)
8. MINIMAP2_SPLICE — all-on-all spliced alignment (source features × target genomes)
9. BAM_TO_GFF — convert spliced BAM alignments to GFF3 gene models
10. GFFCOMPARE — compare mapped vs reference (when target is also source)
11. GFFREAD_EXTRACT — extract mapped transcript sequences
12. BUSCO — completeness assessment (both source and mapped transcripts)
13. REPORT — MultiQC with gffcompare + BUSCO results

### Key Parameters

- `--feature_types` (default: `lnc_RNA,mRNA`) — feature types to extract
- `--minimap2_preset` (default: `splice`) — minimap2 -x preset
- `--busco_lineage` (required) — BUSCO lineage dataset
- `--skip_decoy`, `--skip_busco`, `--skip_gffcompare`, `--skip_download` — skip options

### Reference Scripts

The `scripts/` directory contains the original shell scripts that this pipeline formalizes:
- `scripts/run_minimap_base.sh` — minimap2 splice alignment template
- `scripts/get_allonall_commands.sh` — all-on-all mapping command generation
- `scripts/get_non_overlapping_genes.sh` — bedtools + AGAT non-overlap filtering
- `scripts/get_intergenic_intervals.sh` — AGAT intergenic region extraction
- `scripts/get_relocate_commands.sh` — decoy lncRNA relocation
- `scripts/get_gffcompare_commands.sh` — gffcompare comparison logic
