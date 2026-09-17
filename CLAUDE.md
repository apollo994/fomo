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

**Self-pairs (X→X, only possible with `both`) are run, kept in `allModels`, and excluded
from the top-N ranking pool only.** The projection and its gffcompare are a useful Sn/Pr
ceiling control. Since the single-mapping redesign (plans/15) every source — self included —
is merged into one FASTA and aligned in one run, so self is part of `allModels.raw` and of
its collapse by construction; **`allModels` means literally every source**. The one place it
is dropped is the ranking pool in `consensus_top.nf`: a species' own annotation projected
onto itself is a near-perfect copy, so it would always take the #1 F1 slot and `top3` would
become a restatement of the target's existing annotation. That exclusion lives in exactly one
filter (`meta.id != meta.target_id`, `consensus_top.nf`) — grep for it before adding another.
Consequence: `allModels` accuracy on an annotated `both` target is pulled up towards the self
ceiling and is **not** comparable to the pre-plans/15 `combined` numbers.

Self-projection is also what surfaced `TD2_PREDICT`'s second guard: a self-projection's
lncRNA are exactly the transcripts the upstream PREPROCESSING TD2 pass already cleared of
complete ORFs, so `TD2.LongOrfs --complete-orfs-only` finds **zero** ORFs, PSAURON writes no
`psauron_score.csv`, and `TD2.Predict` dies in `pandas.read_csv`. `modules/local/td2_predict.nf`
therefore checks `td2_work/longest_orfs.pep` *after* LongOrfs as well as checking the input
FASTA before it — nothing to filter is a normal outcome, not a failure. Since plans/15 the
projected-lncRNA TD2 pass runs **once per target on the pooled all-sources model set**, so a
self-projection no longer arrives alone and the zero-ORF path is rare rather than the norm —
but it still fires when a target's only source is itself, or when nothing projects at all.
Keep both guards. Pooling also means the TD2 call is made over a larger training set than the
old per-source passes, so a handful of coding calls can differ from a pre-plans/15 run.

Results are published per target under `${outdir}/targets/<target>/`, in directories named
for **what** they hold rather than for the process that made it (the rules live in the
"Result publishing" section at the bottom of `conf/modules.config`):

```
targets/<target>/alignment/   <target>.allModels.<gtype>.raw.bam(.csi)
targets/<target>/annotation/  <target>.allModels.<gtype>.raw.gff3         (every source)
                              <target>.allModels.<gtype>.collapsed.gff3
                              <target>.top3.<ft>.raw.gff3                 (annotated targets)
                              <target>.top3.<ft>.collapsed.gff3
                              <target>.from_<source>.<gtype>.projected.gff3
targets/<target>/gffcompare/  *.stats *.tracking *.loci *.tmap *.refmap
targets/<target>/select_top_sources/  <target>.top_sources.csv
targets/<target>/multiqc/     multiqc_report.html + multiqc_report_data/
```

plus **one run-level tree**, a sibling of `targets/`, written by `RUN_SUMMARY`:

```
summary/multiqc/  fomo_run_summary.html + fomo_run_summary_data/
summary/tables/   run_*_mqc.tsv, run_*_mqc.json   (the aggregates as data)
                  run_summary.json                (every number, one file)
```

`<gtype>` is `<feature_type>[.decoy]`. There is **one MultiQC report per target** — source-side
sections are target-agnostic and broadcast into every report — **and one for the run**.

# Pipeline Architecture

Conceptual stages and their **implementation status**:

1. **Collect source annotation** — gather lncRNA/mRNA GFF3 + FASTA from source species. *Deferred* — sources are supplied directly via the samplesheet.
2. **Preprocessing** — prepare spliced sequences and decoy sequences for alignment. ✅ Implemented.
3. **Projection** — one splice-aware minimap2 alignment of the merged all-sources sequences onto each target, split back per source afterwards. ✅ Implemented.
4. **Benchmarking** — gffcompare of projected models vs. target reference annotation. ✅ Implemented.
5. **Validation** — splice-junction validation of projected models. ❌ Not yet implemented (next step, see `BRAINSTORM.md`).
6. **Reporting** — one MultiQC report per target, plus one run-level summary report. ✅ Implemented.

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
`F·D` `MERGE_SOURCE_FASTA`, `S` `GUNZIP_FASTA`, `T` `GUNZIP_TARGET`,
`T` `MINIMAP2_INDEX`, **`F·D·T`** minimap2 alignments / `BAM_TO_GFF` /
`FILTER_ALLMODELS` / `SPLIT_GFF_BY_SOURCE` (each split emitting S files),
`T` projected-lncRNA TD2 passes, `F·D·T` `GFFCOMPARE_COMBINE`,
`(F·S·D + 2·F·D)·T_g` benchmark `GFFCOMPARE` **comparisons**, `2·F·T_g` `GFFCOMPARE_TOP`s,
`F·T_g` target references, `T` MultiQC reports, and — independent of every one of those
letters — **one** `RUN_SUMMARY_TABLES` + **one** `MULTIQC_RUN_SUMMARY` per run.

The alignment count is `F·D·T` and **not** `F·S·D·T`: since plans/15 every source's spliced
transcripts are merged into one FASTA per gene type and aligned in a single minimap2 run per
target, then split back apart on the `source=` attribute. The `2·` on the benchmark counts is
the raw/collapsed pair for each aggregate.

**Since plans/18, `GFF_STATS_PROJECTED_BATCH` and `GFFCOMPARE_BATCH` run `F·D·T` and
`F·D·T_g` *tasks* respectively — not `S·T` and `S·T_g`.** The `S·T`/`(F·S·D+2·F·D)·T_g`
figures above (and the `F_S` per-source file counts elsewhere in this doc) still describe
how many per-source *outputs* exist — every source still gets its own stats row and its own
gffcompare `.stats`/`.tracking`/etc. — but one task now produces all of a target's outputs
for that stage in a single SLURM job, looping over every source internally with
`xargs -P ${task.cpus}`. This is the fix for the Nextflow-head OOMs the old S·T/S·T_g task
counts caused on an all-vs-all run of any real size; see `plans/18_per_target_batching.md`.

**There is no gunzip on the GFF3 side at all, and that is deliberate** (plans/17). Every
GFF3 consumer reads `.gz`: `FILTER_TRANSCRIPT`/`FILTER_TARGET` and `GFF_TO_GENE_BED`
decompress inline, and `gff-feature-stats` reads `.gff3.gz` natively — so the samplesheet
GFF3 goes straight into `BENCHMARKING` as well as `PREPROCESSING`. The two surviving
`MAYBE_GUNZIP` aliases are **FASTA-only and exist solely for gffread**, whose `-g` genome
reader has no gzip or bgzf path. `MINIMAP2_INDEX` is fed the gzipped assembly directly
(minimap2 gzopen's its input), so `GUNZIP_TARGET` is off the critical path and feeds only
`EXTRACT_PROJECTED_LNC`. **Do not swap either one for a `samtools faidx` shortcut**: faidx
accepts bgzip but hard-errors on plain gzip, and while the committed test FASTAs are bgzip
(`bin/subsample_test_data.sh`), the real Ensembl/NCBI assemblies are plain gzip — so such a
change passes `-profile test` and fails on the cluster.

Both flags on reproduces the pre-flag behaviour exactly. **`-profile test` sets neither** —
it runs the lncRNA deliverable only (F = 1, D = 1), because both flags together quadruple
the task count and neither is needed to prove the wiring. Pass `--include_mrna` /
`--include_decoy` on the command line to exercise those branches.

Note the asymmetry that makes all-vs-all cheap: **preprocessing scales with S, not S·T**
— it is entirely target-independent (decoys relocate into the *source's own* intergenic
space), so it must never be fanned out per target. `MERGE_SOURCE_FASTA` is target-independent
for the same reason (`F·D` tasks, not `F·D·T`); the fan-out over targets happens *after* it.

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
| `PROJECTION` | **merge every source's spliced FASTA per gene type** (`F·D`, target-independent), index each target, **one minimap2 per (target, gene type)**, BAM→GFF, one projected-lncRNA TD2 pass per target, then **split back per source on `source=`**; one `GFF_STATS_PROJECTED_BATCH` task per target computes every source's stats + both aggregates' (plans/18); gffcompare-collapse `allModels` into `allModels_collapsed` | `gff3` (per target: **List** of every source's projected GFF3 ⊎ both aggregates — plans/18), `allmodels_raw`, `allmodels_collapsed`, `bam`, `index`, `mqc_files` |
| `BENCHMARKING` | filter target ref, **one `GFFCOMPARE_BATCH` task per target** loops gffcompare over every projected model vs. the reference (plans/18), then re-keys each result back to its own source/aggregate; GFF stats on target GFFs. **Targets without a `gff3` are filtered out here** | `stats` (per-source, re-keyed — same shape as before batching), `target_refs`, `mqc_files` |
| `CONSENSUS_TOP` | rank sources by lncRNA transcript F1 **per target**, **subset `allModels` to those sources** (`top3_raw`), gffcompare-collapse it (`top3_collapsed`), score **both** | `gff3`, `stats`, `mqc_files` |
| `REPORTING` | per-target accuracy scatter + one MULTIQC per target | `report`, `data` |
| `RUN_SUMMARY` | digest the pipeline-wide `mqc_files` union + `top_sources` CSVs + a samplesheet-derived roles CSV into run-level tables/plots, then **one MULTIQC for the whole run**. One task per run, no target dimension | `report`, `tables` |

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

**`meta`-map contract:** subworkflows progressively enrich the meta map — `feature_type` + `decoy` (PREPROCESSING), `gtype` = `feature_type[_decoy]` + `target_id` (PROJECTION; `target_id` is also stamped on target-side stats in BENCHMARKING, where it equals `meta.id`), `kind` ∈ {raw, filtered, decoy, projected} (stat producers). **`id` carries a pseudo-source between the merge and the split**: it is the literal `'allModels'` from `MINIMAP2_ALIGN` through `FILTER_ALLMODELS`, becomes the real species again when `SPLIT_GFF_BY_SOURCE`'s output is re-keyed, and is one of `allModels_raw | allModels_collapsed | top3_raw | top3_collapsed` on the four aggregate models. Those four ids are excluded from the ranking pool in **three** places that must agree: the channel filter in `consensus_top.nf`, the re-key closure in `benchmarking.nf` that recovers per-source identity from `GFFCOMPARE_BATCH`'s one-shared-meta batch output (plans/18 — `stats.name[pre.size()..-(suf.size()+1)]`, mirroring `projection.nf`'s identical split-file re-keying), and `AGGREGATE_IDS` in **`bin/fomo_stats.py`** — the shared module that also owns the `<target>.from_<source>.<ft>[.decoy]` stats-filename grammar and the Sn/Pr parser, imported by `select_top_sources.py`, `gffcompare_accuracy_mqc.py` and `run_summary_tables.py`. Nextflow puts `bin/` on PATH but not on Python's import path, so each of those does `sys.path.insert(0, Path(__file__).resolve().parent)` first — the directory is staged/mounted whole, so the sibling import works on every executor (verified under Singularity on the cluster). `target_id` is load-bearing twice over: it is the per-target publishDir segment in `conf/modules.config`, and REPORTING routes MultiQC files by its *presence* (`meta.target_id != null` → that target's report only; absent → broadcast to every report). Downstream modules and `ext.prefix`/`ext.sample_name` closures depend on these keys; preserve them when adding wiring. `feature_type` ranges over the *enabled* subset (see **Feature tracks**), and `decoy: true` / `kind: 'decoy'` occur only with `--include_decoy` — every `meta.decoy` dereference in `conf/modules.config` is Groovy-null-safe and the real track always carries `decoy: false`, so the `ext.prefix`/`ext.sample_name` closures need no change when a track is off. Their `1_raw / 2_lncRNA / 3_decoy_lncRNA / 4_mRNA / 5_decoy_mRNA` ordering ladders are deliberately sparse-safe: a disabled class simply produces no row, and the `'unknown'` fallback cannot trigger because the values are always a subset.

**Filters:** `FILTER_TRANSCRIPT` (always on) keeps only multi-exon (spliced) transcripts. `RELOCATE_LOCI` runs **only with `--include_decoy`**, and caps decoys at `params.decoy_cap` (default 1000; 0 disables) using `params.relocate_seed` for reproducibility.

### The merge/split contract

Merging every source into one alignment is lossless because the species is already carried
end to end, and **nothing in the chain may break that chain of custody**:

`modules/local/rename_fasta_headers.nf` writes `>tid|type|species` → minimap2 puts that in
the BAM QNAME → `bin/bam_to_gff.sh` stamps it on every `gene` and `transcript` row as
`source=<species>` → `bin/gff_by_source.py` keys on that attribute (child `exon`/CDS rows
inherit it through `Parent`). One script serves both directions: split mode
(`--name-template` with a `{source}` placeholder, one file per source) and subset mode
(`--keep <csv> --name`, used for `top3_raw`).

Four things that will bite:

- **The split filename template and the re-keying arithmetic must agree.** `ext.args` for
  `SPLIT_GFF_BY_SOURCE` builds `<target>.from_{source}.<gtype>.projected.gff3`. Since
  plans/18 that per-source list is never re-keyed in `projection.nf` itself (it's folded,
  unexploded, into `GFF_STATS_PROJECTED_BATCH`'s input list) — the re-keying now happens one
  step further downstream, in `benchmarking.nf`, which strips the analogous
  `<target>.from_` / `.<gtype>[.decoy].gffcompare.stats` prefix/suffix *by length* off
  `GFFCOMPARE_BATCH`'s output filenames (not by regex — species names are full of `_` and
  digits) to recover each source's identity for `CONSENSUS_TOP`'s ranking. Change the split
  template and that re-key breaks too, with an out-of-range substring — it does not silently
  mismatch, which is the point.
- **Transcript IDs must be unique across the merged file.** Per-source GFF3s were separate
  files, so a shared id never collided; in one file a duplicate `ID=` makes `gffread -w` and
  the gffcompare collapse silently misbehave. `gff_by_source.py` asserts uniqueness and
  exits 1 naming the id and both sources. Ensembl ids are species-scoped so it should never
  fire.
- **`transpose()` handles both output shapes.** The `path("*.gff3")` glob yields a List in
  split mode and a bare Path in subset mode; `transpose()` passes a non-List element through
  unchanged (verified), so no `arity` declaration is needed.
- **A target with nothing projected onto it is real, not broken — and `SPLIT_GFF_BY_SOURCE`'s
  output is `optional: true` because of it.** If `allModels.raw.gff3` has zero records
  (nothing from any source survived alignment + TD2 filtering onto this target — seen for
  real on divergent targets in a large all-vs-all run), `gff_by_source.py` writes zero files
  and exits 0 by design. `COMBINED_GTF_TO_GFF` then never runs either, since gffcompare's own
  combine-mode output is already `optional: true` upstream. `projection.nf` handles both with
  `join(..., remainder: true)` rather than a plain inner join, which would otherwise make the
  target silently vanish from `GFF_STATS_PROJECTED_BATCH` and benchmarking instead of
  surfacing an honest "0 models" result.

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
  `.gff3.gz` directly, so no gunzip step is needed upstream — on the source side
  (`PREPROCESSING`) or the target side (`BENCHMARKING`).
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
    `'.*:GFF_STATS(_TARGET|_PROJECTED(_BATCH)?|_TOP)?$'` — keep it in sync if a new
    `GFF_STATS` alias is added. **`GFF_STATS_PROJECTED_BATCH`** (plans/18) is the one
    exception to the `ubuntu:22.04` base: it also needs python3 (to run
    `gff_stats_to_mqc.py` in the same task, see Reporting conventions), so it runs on
    `python:3.11` instead. **Not yet verified** to meet the glibc ≥ 2.34 floor — check
    with `docker run --rm python:3.11 ldd --version` before relying on it.
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

The top-level workflow mixes `mqc_files` across subworkflows into `ch_all_mqc` and passes
that one union to **two** consumers — `REPORTING` (one report per target) and
`RUN_SUMMARY` (one report for the run). A channel read twice is forked by Nextflow, so
neither sees a shortened stream. Those are the only two places that call `MULTIQC`
(pinned to **v1.35**; bump in `modules/nf-core/multiqc/{main.nf,environment.yml}` and in
`modules/local/{samtools_to_mqc,select_top_sources}.nf`, which reuse the MultiQC
container. `GFF_STATS_TO_MQC` and `RUN_SUMMARY_TABLES` do not — they need only stdlib
`json`/`csv`, so they run on `python:3.11`).

**The run-level report (`RUN_SUMMARY`)** is deliberately **aggregate-only**: one
`RUN_SUMMARY_TABLES` task digests the union into ~12 run-level tables/plots, and only
those reach its MULTIQC. Do not "simplify" it by handing the raw union to a second
MULTIQC — at S=T=20 that is ~1000 files / ~900 rows (a *bigger* report than the
per-target ones), and every per-target accuracy scatter carries the same hardcoded
section id `gffcompare_accuracy_custom`, so the T of them would overwrite each other.
Two further things it depends on:
- `RUN_SUMMARY` is a **top-level** subworkflow with an **aliased** MULTIQC
  (`MULTIQC_RUN_SUMMARY`). Both matter: `withName` matching is a regex *find*, so
  `FOMO:REPORTING:MULTIQC` (now `$`-anchored, defensively) would otherwise claim a
  same-subworkflow alias, and the `.*:MULTIQC$` publishDir rule is `$`-anchored so the
  alias escapes it — its meta has no `target_id` and would publish to
  `targets/null/multiqc`. Its own publishDir paths are **static strings** for the same
  reason.
- Its plots are **self-describing `*_mqc.json`** (id/section_name/plot_type/pconfig
  embedded, no `sp` pattern, no `custom_data` block); only the TSV tables are declared in
  `assets/multiqc/run_summary_sections.yml`. MultiQC 1.35 accepts a **multi-dataset
  bargraph** (`pconfig.data_labels`) but **not** a multi-dataset heatmap — a list of
  matrices fails validation and the whole report comes out as "No analysis results
  found", so the source × target accuracy matrix is three single-matrix sections
  (F1/Sn/Pr), not one with a switcher.

**One report per target.** `REPORTING` splits the union on the *presence* of
`meta.target_id`: set → the file belongs to that target's report alone; absent → it is
source-side and target-agnostic (input GFF stats, SeqKit, input TD2) and is broadcast into
every report. So a new stat producer only needs to carry `target_id` — or not — and it
routes itself.

**Alignment metrics are per (target, gene type), not per source.** There is one merged BAM,
so `SAMTOOLS_STATS` / `SAMTOOLS_TO_MQC` emit one row per gene type per target
(`<target>.allModels.<1_lncRNA|…>`) — the per-source alignment rows of the pre-plans/15
report are gone, and that is accepted, not a bug. Per-source *model* counts are unaffected:
`GFF_STATS_PROJECTED_BATCH` still computes stats for every split GFF3 — one task per
target now, looping over every source internally, rather than one task per source
(plans/18) — and per-source accuracy still comes from the per-source rows
`GFFCOMPARE_BATCH` re-keys out of its own per-target batch. Recovering per-source
alignment metrics would mean splitting the BAM, which the design deliberately avoids.

Two mechanical traps in that fan-out, both hit for real:
- **Flatten before you `combine`.** An adapter may emit a *list* of paths per item
  (`GFF_STATS_TO_MQC` emits a transcript table *and* a gene-category table), and `combine`
  SPREADS a List-valued item across the output tuple. `[meta, [a, b]].combine(ids)` becomes
  `[a, b, tid]` and the next closure is called with three arguments. `flatMap` to one file
  per item first. Same reason you cannot `.collect()` the shared files and combine *that*.
- **`meta.target_id` is also the publishDir segment** (`conf/modules.config`), so a
  target-scoped process that loses the key publishes to `targets/null/`.

MultiQC config lives in two files under `assets/multiqc/` (wired via
`params.multiqc_main_config` / `params.multiqc_sections_config`, passed as a 2-item
list to MULTIQC's single config slot):
- `main.yml` — top-level layout: `report_title`, `report_comment`, `extra_fn_clean_exts`,
  `report_section_order`, `skip_generalstats: true`, `remove_sections`.
- `sections.yml` — `custom_data` blocks + `sp` patterns for every custom section.

The run-level report has its **own** pair, same split, wired via
`params.multiqc_run_summary_main_config` / `params.multiqc_run_summary_sections_config`:
`run_summary_main.yml` + `run_summary_sections.yml`. The two pairs never mix (separate
tasks, disjoint `run_*_mqc.*` filename suffixes).

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
