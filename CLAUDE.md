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

Conceptual stages and their **implementation status**:

1. **Collect source annotation** — gather lncRNA/mRNA GFF3 + FASTA from source species. *Deferred* — sources are supplied directly via the samplesheet.
2. **Preprocessing** — prepare spliced sequences and decoy sequences for alignment. ✅ Implemented.
3. **Projection** — minimap2 splice-aware alignment of source sequences onto target. ✅ Implemented.
4. **Benchmarking** — gffcompare of projected models vs. target reference annotation. ✅ Implemented.
5. **Validation** — splice-junction validation of projected models. ❌ Not yet implemented (next step, see `BRAINSTORM.md`).
6. **Reporting** — MultiQC report. ✅ Implemented.

### Subworkflow map

The top-level `workflows/fomo.nf` wires four subworkflows under `subworkflows/local/`:

| Subworkflow | Does | Emits |
|-------------|------|-------|
| `PREPROCESSING` | spliced-only filter, spliced-FASTA extraction, intergenic BED, decoy relocation, header rename; GFF + SeqKit stats | `spliced_fasta`, `decoy_spliced_fasta`, `filtered_gff3`, `decoy_gff3`, `intergenic_bed`, `mqc_files` |
| `PROJECTION` | minimap2 align (source spliced FASTA → target), BAM→GFF; samtools stats + GFF stats on projected models | `gff3`, `bam`, `index`, `combined_gff`, `mqc_files` |
| `BENCHMARKING` | filter target ref, gffcompare projected vs. reference; GFF stats on target GFFs | `stats`, `target_refs`, `mqc_files` |
| `REPORTING` | the single MULTIQC invocation | `report`, `data` |

**Single-target assumption:** `PROJECTION` calls `.first()` on the index channel, so it indexes/aligns against exactly one target. The samplesheet schema permits more, but only the first target is used. Generalise here if multi-target support is needed.

**`meta`-map contract:** subworkflows progressively enrich the meta map — `feature_type` + `decoy` (PREPROCESSING), `target_id` (PROJECTION), `kind` ∈ {raw, filtered, decoy, projected} (stat producers). Downstream modules and `ext.prefix`/`ext.sample_name` closures depend on these keys; preserve them when adding wiring.

**Always-on filters:** `FILTER_TRANSCRIPT` keeps only multi-exon (spliced) transcripts; `RELOCATE_LOCI` caps decoys at `params.decoy_cap` (default 1000; 0 disables) using `params.relocate_seed` for reproducibility.

### Preprocessing detail

Each step maps to a legacy script in `legacy_scripts/preprocessing/` which serves as the reference implementation:

| Step | Legacy script | Tool | Purpose |
|------|--------------|------|---------|
| Extract lncRNA/mRNA spliced FASTA | `00_get_feature_annotation.sh` | AGAT | Filter GFF3 by feature type, keep longest isoforms, extract exon sequences |
| Extract intergenic intervals | `02_get_intergenic_intervals.sh` | AGAT + bedtools | Build intergenic BED from target annotation |
| Build decoy sequences | `03_relocate_loci.py` | Python | Relocate source loci into intergenic intervals of target to create decoy FASTA |
| Extract decoy spliced FASTA | `04_get_decoy_sequence_commands.sh` | AGAT | Extract exon sequences from decoy GFF3 |
| Preprocessing statistics | `05_get_gff_statistics_commands.sh` | AGAT (now `gff-feature-stats`) | Collect GFF stats at each stage |

The projection stage uses minimap2 (see `legacy_scripts/minimap_transfer/`) and converts BAM → GFF3 with exon structure.

# Key Tools
- **gffread** - GFF3 manipulation. Prefer this over **AGAT** when possible. 
- **gff-feature-stats** — GFF3 feature statistics (gene categories, per-transcript-type
  counts/lengths, introns). Replaced `agat_sp_statistics.pl` for every stats step: same
  numbers (verified metric-by-metric against AGAT, introns included), plus intron stats
  AGAT-style parsing never gave us, at a flat ~20 MB RSS instead of 8–12 GB. Reads
  `.gff3.gz` directly, so no gunzip step is needed upstream.
  **Vendored as `bin/gff-feature-stats`** — a release build of v0.2.0 from
  `github.com/apollo994/gff-feature-stats` @ `67ffea7`. The **commit is the provenance**:
  `67ffea7` changed behaviour without bumping the crate version, so `-V` (and the
  versions topic) reports `0.2.0` for it and the earlier `6b31e3a` build alike. It is
  dynamically linked and needs glibc ≥ 2.34, which the pinned `ubuntu:22.04` task image
  satisfies; re-test before pointing `GFF_STATS` at an older base image. To rebuild (bump
  the commit in `modules/local/gff_stats.nf` when you do):
  ```sh
  cd ~/repos/gff-feature-stats && cargo build --release   # or, for a static binary:
  RUSTFLAGS='-C target-feature=+crt-static' cargo build --release
  cp target/release/gff-feature-stats <fomo>/bin/ && chmod 755 <fomo>/bin/gff-feature-stats
  singularity exec docker://ubuntu:22.04 ./bin/gff-feature-stats -V   # smoke test
  ```
  Three constraints to know:
  - **Input must be grouped by seqid** (contiguous records per seqid; order *within* a
    seqid is free). Interleaved seqids exit 1, but only above the tool's 200k-line batch
    threshold — smaller files form one batch and are always accepted, so a small test run
    can pass where a full-size one fails. `bin/relocate_loci.py` sorts its decoys by
    seqid for exactly this reason; every other GFF here is already grouped (Ensembl
    input, gffread output, and `BAM_TO_GFF` output via `samtools sort`).
  - It cannot compute genome coverage (no genome-size input).
  - It resolves a transcript's gene via `Parent`/`gene`/`Gene` only — so gffread-derived
    GFF3s (`combined`, `top3`), which carry the gene in `geneID=` and emit no gene
    features, report 0 genes.
- **AGAT** — GFF3 manipulation (filtering, longest isoform selection). No longer used by the
  pipeline; kept as the reference implementation in `legacy_scripts/`.
- **minimap2** — splice-aware long-read alignment (source spliced FASTA → target)
- **bedtools** — genomic interval arithmetic
- **samtools** — BAM handling
- **gffcompare** — annotation comparison for benchmarking

# Reporting conventions

Every subworkflow that emits statistics:
- Runs its stat producers (`gff-feature-stats`, SeqKit, samtools stats, gffcompare, ...)
  and any required adapter modules (modules named `*_TO_MQC` that parse a tool's output
  into an `*_mqc.tsv` — usually one row, but `GFF_STATS_TO_MQC` emits one row per
  transcript type and a second file for gene categories).
- Emits a `mqc_files` channel of shape `tuple(meta, path)` — same shape as the rest
  of the pipeline so the meta is available for tracing.
- Subworkflows with no stats emit `Channel.empty()` as `mqc_files`.

The top-level workflow mixes `mqc_files` across subworkflows and passes the union to
`REPORTING`. `REPORTING` strips meta, collects paths, and is the only place that calls
`MULTIQC` (pinned to **v1.35**; bump in `modules/nf-core/multiqc/{main.nf,environment.yml}`
and in `modules/local/{samtools_to_mqc,select_top_sources}.nf`, which reuse the MultiQC
container. `GFF_STATS_TO_MQC` does not — it needs only stdlib `json`, so it runs on
`python:3.11`).

MultiQC config lives in two files under `assets/multiqc/` (wired via
`params.multiqc_main_config` / `params.multiqc_sections_config`, passed as a 2-item
list to MULTIQC's single config slot):
- `main.yml` — top-level layout: `report_title`, `report_comment`, `extra_fn_clean_exts`,
  `report_section_order`, `skip_generalstats: true`, `remove_sections`.
- `sections.yml` — `custom_data` blocks + `sp` patterns for every custom section.

### Adding a new custom-content section — rules learned the hard way

1. **Route by `sp` pattern only, never both `sp` and an embedded `# id:` header.**
   Each adapter emits a file with a **unique suffix** (`*_gffstats_input_mqc.tsv`,
   `*_gffstats_genes_mqc.tsv`, `*_gffstats_projection_mqc.tsv`, `*_seqkit_mqc.tsv`,
   `*_samtools_align_mqc.tsv`) matched
   1:1 by an `sp.<section>.fn` pattern. The `*_TO_MQC` adapters deliberately do **not**
   write a `# id:` header — a file discovered by both `sp` *and* a header renders the
   section **twice**.
2. **Order sections with `report_section_order` (numeric `order`), not
   `custom_content.order`.** Combining `custom_content.order` with default discovery
   also double-renders every listed section.
3. **Native modules' tables are not configurable column-by-column.** To show a curated
   column set (e.g. the samtools alignment metrics), write a small `*_TO_MQC` adapter
   that emits exactly the wanted columns as a custom-content table, then hide the native
   section via `remove_sections` (samtools' `samtools-stats` violin is removed this way;
   its "Percent mapped" bar chart is kept).

When a new tool produces statistics: add a `custom_data` block + a uniquely-suffixed
`sp` pattern to `sections.yml`, give it a slot in `report_section_order` in `main.yml`,
and add any new filename suffixes to `extra_fn_clean_exts`.

# Parameters & validation

Pipeline params live in `nextflow.config` and are validated against `nextflow_schema.json`
by `validateParameters()` (nf-schema) at the top of `workflows/fomo.nf`. Validation runs
with `failUnrecognisedParams = true`.

**Therefore: any new param added to `nextflow.config` MUST also be added to
`nextflow_schema.json`, or the pipeline aborts at launch.** Keep the two in sync.
Profile-metadata params that aren't pipeline inputs go in
`validation.defaultIgnoreParams` instead of the schema. `nextflow run . --help` renders
the schema as a grouped help menu. The samplesheet (not params) is validated separately
by `assets/schema_input.json`.
