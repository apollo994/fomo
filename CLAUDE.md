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

# Samplesheet & roles

Columns: `species,role,fasta,gff3` — validated by `assets/schema_input.json`.

`role` ∈ `source` | `target` | `both`. **`both` is the all-vs-all case**: the species
donates its annotation *and* receives projections, so one run replaces the N runs with N
hand-written samplesheets that the pre-multitarget pipeline needed. A species is one row —
`species` must be unique, and duplicates abort at launch (it keys every join).

`gff3` is **required for `source` and `both`** (enforced by an `allOf`/`if`-`then` in the
schema) and **optional for a pure `target`**: an un-annotated assembly is projected onto
and gets per-source models plus the all-source `combined` consensus, but no gffcompare, no
target GFF stats, and no top-N consensus (the ranking is derived from gffcompare F1). This
falls out structurally — `BENCHMARKING` filters to targets with a GFF and its projection ⋈
reference join is an inner join, so nothing downstream needs a guard.

**Self-pairs (X→X, only possible with `both`) are run but excluded from consensus.** The
projection and its gffcompare are a useful Sn/Pr ceiling control, but a species' own
annotation projected onto itself is a near-perfect copy, so it is dropped from the
all-source `combined` model (`projection.nf`) and from the top-N ranking pool
(`consensus_top.nf`) — otherwise it would always rank #1 on F1 and the consensus would be a
restatement of the target's existing annotation. Consequence: a target whose only source is
itself produces no `combined` model.

Self-projection also drives `TD2_PREDICT`'s second guard: its projected lncRNA are exactly
the transcripts the upstream PREPROCESSING TD2 pass already cleared of complete ORFs, so
`TD2.LongOrfs --complete-orfs-only` finds **zero** ORFs, PSAURON writes no
`psauron_score.csv`, and `TD2.Predict` dies in `pandas.read_csv`. `modules/local/td2_predict.nf`
therefore checks `td2_work/longest_orfs.pep` *after* LongOrfs as well as checking the input
FASTA before it — nothing to filter is a normal outcome, not a failure.

Results are published per target under `${outdir}/targets/<target>/<process>/`, and there
is **one MultiQC report per target** (`targets/<target>/multiqc/`) — source-side sections
are target-agnostic and broadcast into every report.

# Pipeline Architecture

Conceptual stages and their **implementation status**:

1. **Collect source annotation** — gather lncRNA/mRNA GFF3 + FASTA from source species. *Deferred* — sources are supplied directly via the samplesheet.
2. **Preprocessing** — prepare spliced sequences and decoy sequences for alignment. ✅ Implemented.
3. **Projection** — minimap2 splice-aware alignment of source sequences onto target. ✅ Implemented.
4. **Benchmarking** — gffcompare of projected models vs. target reference annotation. ✅ Implemented.
5. **Validation** — splice-junction validation of projected models. ❌ Not yet implemented (next step, see `BRAINSTORM.md`).
6. **Reporting** — MultiQC report. ✅ Implemented.

### Feature tracks (what gets transferred)

Only **lncRNA** is transferred by default — it is the deliverable. Two development
instruments are opt-in:

| Flag | Adds | Why |
|------|------|-----|
| `--include_mrna` | the source `mRNA` track through every stage | **positive control** — mRNA projects well from close relatives, so its gffcompare Sn/Pr bound what the lncRNA track can reach |
| `--include_decoy` | a relocated decoy track per *enabled* feature type | **negative control / FDR baseline** — source loci moved into the source's own intergenic space, so anything that still projects is a false positive |

With **S** sources (`role` ∈ {source, both}), **T** targets (`role` ∈ {target, both}),
**T_g** ≤ T of them carrying a `gff3`, F enabled feature types (1, or 2 with
`--include_mrna`) and D = 2 with `--include_decoy` else 1: `F·S` `FILTER_TRANSCRIPT`,
`T` `MINIMAP2_INDEX`, `F·S·T·D` minimap2 alignments, `F·T·D` `GFFCOMPARE_COMBINE` groups,
`(F·S·D + F·D)·T_g` benchmark `GFFCOMPARE`s, `F·T_g` target references, `T` MultiQC reports.
Both flags on reproduces the pre-flag behaviour exactly. **`-profile test` sets neither** —
it runs the lncRNA deliverable only (F = 1, D = 1), because both flags together quadruple
the task count and neither is needed to prove the wiring. Pass `--include_mrna` /
`--include_decoy` on the command line to exercise those branches.

Note the asymmetry that makes all-vs-all cheap: **preprocessing scales with S, not S·T**
— it is entirely target-independent (decoys relocate into the *source's own* intergenic
space), so it must never be fanned out per target.

The enabled list is derived **once**, in `workflows/fomo.nf`
(`ch_feature_types = Channel.value(['lnc_RNA'] + (params.include_mrna ? ['mRNA'] : []))`),
and passed as a `take:` input to both `PREPROCESSING` (source filtering) and `BENCHMARKING`
(target-reference filtering). **Never re-derive it inside a subworkflow** — benchmarking
pairs projections to references with an inner `combine(by: 0)` on `feature_type`, so a list
that disagrees with the source side silently *drops* projections instead of failing.

**Decoy placement is annotation-wide, not track-wide.** `GFF_TO_GENE_BED` reads the *raw*
samplesheet GFF3 (not the per-feature-type `FILTER_TRANSCRIPT` output) and keeps every
gene-level feature (`gene`, `pseudogene`, `ncRNA_gene`), so `BEDTOOLS_COMPLEMENT` yields the
complement of the **entire** annotation — coding genes included — regardless of
`--include_mrna`. The chain runs once per source (N tasks, not `F·N`) and its single BED is
broadcast to every feature type. So: gate it on `include_decoy` only, and never repoint
`GFF_TO_GENE_BED` at a filtered GFF3, which would redefine intergenic as "outside lncRNA
genes" and let lncRNA decoys land inside real mRNA loci.

### Subworkflow map

The top-level `workflows/fomo.nf` wires five subworkflows under `subworkflows/local/`.
It also owns the role split: `role` ∈ {source, both} → `ch_sources`, `role` ∈ {target,
both} → `ch_targets`. Two `.filter`s, **not** `.branch` — branch routes each item to
exactly one output, so a `both` row could never reach both channels.

| Subworkflow | Does | Emits |
|-------------|------|-------|
| `PREPROCESSING` | spliced-only filter, spliced-FASTA extraction, intergenic BED, decoy relocation, header rename; GFF + SeqKit stats. **Per source, target-independent — S tasks, not S·T** | `spliced_fasta`, `decoy_spliced_fasta`*, `filtered_gff3`, `decoy_gff3`*, `intergenic_bed`*, `mqc_files` |
| `PROJECTION` | index each target, minimap2 align every source track × every target, BAM→GFF; samtools stats + GFF stats on projected models; per-(target, gene type) `combined` consensus | `gff3`, `bam`, `index`, `combined_gff`, `mqc_files` |
| `BENCHMARKING` | filter target ref, gffcompare projected vs. reference; GFF stats on target GFFs. **Targets without a `gff3` are filtered out here** | `stats`, `target_refs`, `mqc_files` |
| `CONSENSUS_TOP` | rank sources by lncRNA transcript F1 **per target**, gffcompare-combine the top N, score the result | `gff3`, `stats`, `mqc_files` |
| `REPORTING` | per-target accuracy scatter + one MULTIQC per target | `report`, `data` |

\* `Channel.empty()` without `--include_decoy` — the whole intergenic/relocation branch
(`SAMTOOLS_FAIDX` → `GFF_TO_GENE_BED` → `BEDTOOLS_COMPLEMENT` → `RELOCATE_LOCI` →
`EXTRACT_DECOY_SEQUENCES`) is skipped by an `if (params.include_decoy)` guard in
`preprocessing.nf`.

**Multi-target wiring — two rules.**

*1. Never nest the subworkflows.* The obvious way to add targets is a `PER_TARGET`
subworkflow looped over targets. It breaks ~30 `ext.prefix` / `ext.sample_name` closures
**silently**: every selector in `conf/modules.config` is a fully-qualified process path
(`'FOMO:PROJECTION:MINIMAP2_ALIGN'`, …) and `withName` matching is a regex *find*, so an
extra level makes the substring stop matching and modules fall back to their in-module
defaults — `${meta.id}`, i.e. the **source** id, producing exactly the collisions those
closures exist to prevent. The target dimension lives in the *channels*, not the call tree.

*2. Every projection ⋈ reference join keys on `[target_id, feature_type]`, never
`feature_type` alone* (`benchmarking.nf`, `consensus_top.nf`). On feature_type alone each
projection is also compared against every *other* target's annotation, and since
`ext.prefix` is built from the **query** meta, all T of those tasks emit the same filename
into the same published directory — wrong numbers under a plausible-looking name. The same
applies to `SELECT_TOP_SOURCES`: it is grouped per target because
`bin/select_top_sources.py` keys its score map by source alone, so pooling targets would
have them overwrite each other.

**`meta`-map contract:** subworkflows progressively enrich the meta map — `feature_type` + `decoy` (PREPROCESSING), `target_id` (PROJECTION; also stamped on target-side stats in BENCHMARKING, where it equals `meta.id`), `kind` ∈ {raw, filtered, decoy, projected} (stat producers). `target_id` is load-bearing twice over: it is the per-target publishDir segment in `nextflow.config`, and REPORTING routes MultiQC files by its *presence* (`meta.target_id != null` → that target's report only; absent → broadcast to every report). Downstream modules and `ext.prefix`/`ext.sample_name` closures depend on these keys; preserve them when adding wiring. `feature_type` ranges over the *enabled* subset (see **Feature tracks**), and `decoy: true` / `kind: 'decoy'` occur only with `--include_decoy` — every `meta.decoy` dereference in `conf/modules.config` is Groovy-null-safe and the real track always carries `decoy: false`, so the `ext.prefix`/`ext.sample_name` closures need no change when a track is off. Their `1_raw / 2_lncRNA / 3_decoy_lncRNA / 4_mRNA / 5_decoy_mRNA` ordering ladders are deliberately sparse-safe: a disabled class simply produces no row, and the `'unknown'` fallback cannot trigger because the values are always a subset.

**Filters:** `FILTER_TRANSCRIPT` (always on) keeps only multi-exon (spliced) transcripts. `RELOCATE_LOCI` runs **only with `--include_decoy`**, and caps decoys at `params.decoy_cap` (default 1000; 0 disables) using `params.relocate_seed` for reproducibility.

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
  - **It is an x86-64 binary, and `ubuntu:22.04` is multi-arch.** On an arm64 host
    (Apple Silicon) Docker pulls the arm64 image, which has no
    `/lib64/ld-linux-x86-64.so.2`, and every `GFF_STATS*` task dies with
    `rosetta error: failed to open elf ...` → exit **133**. The `docker` profile in
    `nextflow.config` therefore pins these processes to `--platform=linux/amd64`
    (no-op on x86-64 hosts; scoped to that profile so the Singularity/HPC path is
    untouched). The selector is
    `'.*:GFF_STATS(_TARGET|_PROJECTED|_TOP)?$'` — keep it in sync if a new
    `GFF_STATS` alias is added.
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
`REPORTING`, which is the only place that calls `MULTIQC` (pinned to **v1.35**; bump in
`modules/nf-core/multiqc/{main.nf,environment.yml}` and in
`modules/local/{samtools_to_mqc,select_top_sources}.nf`, which reuse the MultiQC
container. `GFF_STATS_TO_MQC` does not — it needs only stdlib `json`, so it runs on
`python:3.11`).

**One report per target.** `REPORTING` splits the union on the *presence* of
`meta.target_id`: set → the file belongs to that target's report alone; absent → it is
source-side and target-agnostic (input GFF stats, SeqKit, input TD2) and is broadcast into
every report. So a new stat producer only needs to carry `target_id` — or not — and it
routes itself.

Two mechanical traps in that fan-out, both hit for real:
- **Flatten before you `combine`.** An adapter may emit a *list* of paths per item
  (`GFF_STATS_TO_MQC` emits a transcript table *and* a gene-category table), and `combine`
  SPREADS a List-valued item across the output tuple. `[meta, [a, b]].combine(ids)` becomes
  `[a, b, tid]` and the next closure is called with three arguments. `flatMap` to one file
  per item first. Same reason you cannot `.collect()` the shared files and combine *that*.
- **`meta.target_id` is also the publishDir segment** (`nextflow.config`), so a
  target-scoped process that loses the key publishes to `targets/null/`.

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
