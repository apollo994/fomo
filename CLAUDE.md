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

`gff3` is **required for `source` and `both`** (an `allOf`/`if`-`then` rule in the schema)
and **optional for a pure `target`** (plans/22; plans/20 had forbidden it). The two uses
of a gff3 are different:

- **Benchmarking** needs the species' own *filtered* annotation (see **The
  filtered-annotation contract**), which only a donor has. So the **benchmarked targets
  are exactly the `both` rows**, and a species cannot be benchmarked without also donating.
- **Curation** (plans/21) uses the target's *samplesheet* gff3 as is — to drop overlapping
  models, tag each curated gene with its reference gene (`ref_gene_*`) and build the merged
  annotation. Any target with a gff3 gets it, `target` or `both`.

So a pure target **with** a gff3 is curated against it and gets its raw GFF stats in its
report and in the run summary, but no gffcompare; `workflows/fomo.nf`
logs these targets at launch. A pure target **without** one is curated reference-free
(`ref_location=no_reference`, no merged file). Neither needs a guard downstream: it has no
filtered annotation, and the projection ⋈ reference join is an inner join on it.

**There is no self-projection (plans/23).** A `both` species donates to every *other*
target but is never aligned to its own genome, on any track (lncRNA, mRNA, decoys). The
merged source FASTA is still built once per gene type and shared by every target (plans/15),
so the exclusion happens at alignment: `projection.nf` sets `exclude_self: true` on the align
meta of a `both` target, and the FOMO-local branch of `modules/nf-core/minimap2/align` drops
the query records whose header species (the last `|`-field of `>tid|type|species`) is the
target, into a temporary FASTA, before minimap2 runs. So **`allModels` means every source
except the target itself**, the per-source split never has a `from_<target>` file, and the
accuracy heatmap has no diagonal. Self-projection used to be the Sn/Pr ceiling control and
inflated `allModels` (up to ~90 % of a target's raw models); it contributed nothing to the
curated deliverable, which never counts a target's own projections as support.
`bin/curate_models.py` still drops `source=<target>` models with a warning, for inputs from
pre-23 runs. When a `both` target's ONLY source is itself, its query is empty: minimap2
writes a header-only BAM (verified) and the "nothing projected" path takes over.

Self-projection is also what surfaced `TD2_PREDICT`'s second guard (kept): a self-projection's
lncRNA were exactly the transcripts the upstream PREPROCESSING TD2 pass already cleared of
complete ORFs, so `TD2.LongOrfs --complete-orfs-only` finds **zero** ORFs, PSAURON writes no
`psauron_score.csv`, and `TD2.Predict` dies in `pandas.read_csv`. `modules/local/td2_predict.nf`
therefore checks `td2_work/longest_orfs.pep` *after* LongOrfs as well as checking the input
FASTA before it — nothing to filter is a normal outcome, not a failure. Since plans/15 the
projected-lncRNA TD2 pass runs **once per target on the pooled all-sources model set**, so a
single source's models no longer arrive alone and the zero-ORF path is rare rather than the
norm — but it still fires when nothing (or only ORF-free models) projects at all.
Keep both guards. Pooling also means the TD2 call is made over a larger training set than the
old per-source passes, so a handful of coding calls can differ from a pre-plans/15 run.

Results are published per target under `${outdir}/targets/<target>/`, in directories named
for **what** they hold rather than for the process that made it (the rules live in the
"Result publishing" section at the bottom of `conf/modules.config`):

```
targets/<target>/alignment/   <target>.allModels.<gtype>.raw.bam(.csi)
targets/<target>/annotation/  <target>.allModels.<gtype>.raw.gff3         (every source but the target)
                              <target>.from_<source>.<gtype>.projected.gff3
targets/<target>/gffcompare/  *.stats *.tracking *.loci *.tmap *.refmap   (role 'both' targets)
targets/<target>/curated/     <target>.curated.<gtype>.gff3.gz            (plans/21, every track)
                              <target>.curated.<gtype>.merged.gff3.gz     (annotated targets)
                              <target>.curated.<gtype>.curation.{tsv,json}
targets/<target>/multiqc/     multiqc_report.html + multiqc_report_data/
```

plus **one per-source tree** (plans/20), keyed by the source id, holding each source's one
filtered annotation — spliced, longest isoform per gene, lncRNA TD2-noncoding; exactly the
transcripts that were projected, and the gffcompare reference for a `both` species:

```
sources/<source>/annotation/  <source>.<ft>.filtered.gff3
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
**T_g** ≤ T of them benchmarked (the `both` rows — plans/22), F enabled feature types (1, or 2 with
`--include_mrna`) and D = 2 with `--include_decoy` else 1: `F·S` `FILTER_TRANSCRIPT`,
`F·D` `MERGE_SOURCE_FASTA`, `S` `GUNZIP_FASTA`, `T` `GUNZIP_TARGET`,
`T` `MINIMAP2_INDEX`, **`F·D·T`** minimap2 alignments / `BAM_TO_GFF` /
`FILTER_ALLMODELS` / `SPLIT_GFF_BY_SOURCE` (each split emitting S files, S−1 for a `both`
target — no self-projection), `T` projected-lncRNA TD2 passes,
`(F·(S−1)·D + F·D)·T_g` benchmark `GFFCOMPARE` **comparisons** (every other source + `allModels_raw`),
`F·S` `FILTER_ANNOTATION` (the filtered annotation, which doubles as the reference — there is
no separate target-reference filter since plans/20), `F·D·T` `CURATE_MODELS` (plans/21; none with
`--curate false`), `T` MultiQC reports, and — independent of every one of those
letters — **one** `RUN_SUMMARY_TABLES` + **one** `MULTIQC_RUN_SUMMARY` per run.

The alignment count is `F·D·T` and **not** `F·S·D·T`: since plans/15 every source's spliced
transcripts are merged into one FASTA per gene type and aligned in a single minimap2 run per
target, then split back apart on the `source=` attribute. The `+ F·D` on the benchmark counts
is the one aggregate, `allModels_raw` (the gffcompare-collapsed `allModels_collapsed` and the
top-N consensus were removed in plans/23).

**Since plans/18, `GFF_STATS_PROJECTED_BATCH` and `GFFCOMPARE_BATCH` run `F·D·T` and
`F·D·T_g` *tasks* respectively — not `S·T` and `S·T_g`.** The `S·T`/`(F·(S−1)·D + F·D)·T_g`
figures above (and the `F_S` per-source file counts elsewhere in this doc) still describe
how many per-source *outputs* exist — every source still gets its own stats row and its own
gffcompare `.stats`/`.tracking`/etc. — but one task now produces all of a target's outputs
for that stage in a single SLURM job, looping over every source internally with
`xargs -P ${task.cpus}`. This is the fix for the Nextflow-head OOMs the old S·T/S·T_g task
counts caused on an all-vs-all run of any real size; see `plans/18_per_target_batching.md`.

**There is no gunzip on the GFF3 side at all, and that is deliberate** (plans/17). Every
GFF3 consumer reads `.gz`: `FILTER_TRANSCRIPT` and `GFF_TO_GENE_BED` decompress inline,
and `gff-feature-stats` reads `.gff3.gz` natively — so the samplesheet GFF3 goes straight
into `PREPROCESSING`, and into `BENCHMARKING` for the target's `kind: 'raw'` stats. The two surviving
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
(`feature_types = ['lnc_RNA'] + (params.include_mrna ? ['mRNA'] : [])`, a plain List), and
passed as a `take:` input to `PREPROCESSING` only. `BENCHMARKING` gets its references from
`PREPROCESSING.out.filtered_gff3` — already split by that same list — so its inner
`combine(by: [0, 1])` on `[target_id, feature_type]` cannot disagree with the source side
by construction (plans/20). **Never re-derive the list inside a subworkflow**, and never
give BENCHMARKING its own feature-type filter again: a list that disagrees with the source
side silently *drops* projections instead of failing.

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
| `PREPROCESSING` | spliced-only + longest-isoform filter, spliced-FASTA extraction, header rename (+ id map), input TD2, the one filtered annotation (`FILTER_ANNOTATION`), intergenic BED, decoy relocation; GFF + SeqKit + transcript-filter stats. **Per source, target-independent — S tasks, not S·T** | `spliced_fasta`, `decoy_spliced_fasta`*, `filtered_gff3` (describes `spliced_fasta` exactly; published, and the benchmarking reference), `decoy_gff3`*, `intergenic_bed`*, `mqc_files` |
| `PROJECTION` | **merge every source's spliced FASTA per gene type** (`F·D`, target-independent), index each target, **one minimap2 per (target, gene type)** with the target's own records dropped from the query (`exclude_self`, plans/23), BAM→GFF, one projected-lncRNA TD2 pass per target, then **split back per source on `source=`**; one `GFF_STATS_PROJECTED_BATCH` task per target computes every source's stats + the `allModels_raw` aggregate's (plans/18) | `gff3` (per target: **List** of every source's projected GFF3 ⊎ `allModels_raw` — plans/18), `allmodels_raw`, `bam`, `index`, `mqc_files` |
| `BENCHMARKING` | take each `both` target's reference from `PREPROCESSING.out.filtered_gff3` (no filtering of its own), **one `GFFCOMPARE_BATCH` task per target** loops gffcompare over every projected model vs. the reference (plans/18), then re-keys each result back to its own source/aggregate; GFF stats on target GFFs. **Targets without a `gff3` are filtered out here** | `stats` (per-source, re-keyed — same shape as before batching), `mqc_files` |
| `CURATION` | plans/21. **One `CURATE_MODELS` task per (target, feature type, decoy)** on `allModels.<gtype>.raw.gff3` (`bin/curate_models.py` + `bin/curation_to_mqc.py`, one Wave image with python + gffcompare 0.12.6 + gffread 0.12.7). Keeps intron chains shared **exactly** by ≥ `curate_min_species` species (gffcompare `--no-merge` cliques, cross-checked in Python; self never counts), drops models with an exon on any **same-strand reference exon** (every biotype), links supported chains sharing a splice junction into genes, one representative per gene (most SJs → species → length → de), `ref_location` tag, then merges into the reference. The reference is the target's **samplesheet** GFF3, *not* the filtered ★ annotation; a pure target is curated without one (no merged file). Gated on `params.curate` via `ext.when` on `CURATE_MODELS` (always called, so its selectors never warn) | `curated`, `merged`, `report`, `mqc_files` (3 tables, target-scoped) |
| `REPORTING` | per-target accuracy scatter + one MULTIQC per target | `report`, `data` |
| `RUN_SUMMARY` | digest the pipeline-wide `mqc_files` union (curation tables included) + a samplesheet-derived roles CSV into run-level tables/plots, then **one MULTIQC for the whole run**. One task per run, no target dimension | `report`, `tables` |

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
`feature_type` alone* (`benchmarking.nf`). On feature_type alone each
projection is also compared against every *other* target's annotation, and since
`ext.prefix` is built from the **query** meta, all T of those tasks emit the same filename
into the same published directory — wrong numbers under a plausible-looking name. (The
curation join in `curation.nf` keys on `target_id` alone, on purpose: one reference per
target serves every feature type.)

**`meta`-map contract:** subworkflows progressively enrich the meta map — `feature_type` + `decoy` (PREPROCESSING), `gtype` = `feature_type[_decoy]` + `target_id` (PROJECTION; `target_id` is also stamped on target-side stats in BENCHMARKING, where it equals `meta.id`), `kind` ∈ {raw, filtered, decoy, projected} (stat producers), `exclude_self` (PROJECTION's align meta: the target is also a source, so MINIMAP2_ALIGN drops its own records — plans/23; it rides along downstream harmlessly). **`id` carries a pseudo-source between the merge and the split**: it is the literal `'allModels'` from `MINIMAP2_ALIGN` through `FILTER_ALLMODELS`, becomes the real species again when `SPLIT_GFF_BY_SOURCE`'s output is re-keyed, `'curated'` in CURATION, and `allModels_raw` on the one aggregate model (the batch tasks symlink it under that id). That id is kept out of per-source tables in **two** places that must agree: the re-key closure in `benchmarking.nf` that recovers per-source identity from `GFFCOMPARE_BATCH`'s one-shared-meta batch output (plans/18 — `stats.name[pre.size()..-(suf.size()+1)]`, mirroring `projection.nf`'s identical split-file re-keying), and `AGGREGATE_IDS` in **`bin/fomo_stats.py`** — the shared module that also owns the `<target>.from_<source>.<ft>[.decoy]` stats-filename grammar and the Sn/Pr parser, imported by `gffcompare_accuracy_mqc.py` and `run_summary_tables.py`. Nextflow puts `bin/` on PATH but not on Python's import path, so each of those does `sys.path.insert(0, Path(__file__).resolve().parent)` first — the directory is staged/mounted whole, so the sibling import works on every executor (verified under Singularity on the cluster). `target_id` is load-bearing twice over: it is the per-target publishDir segment in `conf/modules.config`, and REPORTING routes MultiQC files by its *presence* (`meta.target_id != null` → that target's report only; absent → broadcast to every report). Downstream modules and `ext.prefix`/`ext.sample_name` closures depend on these keys; preserve them when adding wiring. `feature_type` ranges over the *enabled* subset (see **Feature tracks**), and `decoy: true` / `kind: 'decoy'` occur only with `--include_decoy` — every `meta.decoy` dereference in `conf/modules.config` is Groovy-null-safe and the real track always carries `decoy: false`, so the `ext.prefix`/`ext.sample_name` closures need no change when a track is off. Their `1_raw / 2_lncRNA / 3_decoy_lncRNA / 4_mRNA / 5_decoy_mRNA` ordering ladders are deliberately sparse-safe: a disabled class simply produces no row, and the `'unknown'` fallback cannot trigger because the values are always a subset.

**Filters:** `FILTER_TRANSCRIPT` (always on) keeps only **stranded** transcripts — own strand `+`/`-`, every child on that same strand, no parent with an undefined strand; this drops RefSeq's trans-spliced models (plant mitochondrial *nad1*/*nad2*/*rps12*, strand `?`), on which `gffread -w` hard-errors ("Error parsing strand") — then only multi-exon (spliced) transcripts and, with `--longest_isoform` (default on, plans/20), only **one per gene: the longest spliced isoform** by summed exon length (ties → smallest transcript ID; selection runs *after* the multi-exon filter, so a gene whose longest isoform is mono-exonic keeps its longest spliced one). It is the **only** place selection happens — the reference is derived from its output, so the two sides cannot drift. It also writes a funnel row (`*_transcript_filter_mqc.tsv`: n_transcripts / n_stranded / n_spliced / n_selected) that feeds the `transcript_filter` MultiQC section and the run summary's `dropped_unstranded` / `dropped_single_exon` / `dropped_non_longest` columns. Decoys are relocated from these **pre-TD2 candidates**, not from the filtered annotation: a decoy borrows only a model's shape and takes its sequence from intergenic DNA, so TD2's verdict on the original locus is irrelevant to it. `RELOCATE_LOCI` runs **only with `--include_decoy`**, and caps decoys at `params.decoy_cap` (default 1000; 0 disables) using `params.relocate_seed` for reproducibility.

### The filtered-annotation contract

Each source has **one** filtered annotation per feature type,
`<source>.<ft>.filtered.gff3` (`FILTER_ANNOTATION`, plans/20), and it is used three ways:
it describes exactly the FASTA sent to projection, it is the gffcompare reference for a
`both` species (and so for the target `kind: 'filtered'` stats), and it
is published under `sources/<source>/annotation/`. It is built **from** the projected FASTA,
not alongside it:

`FILTER_TRANSCRIPT` (candidates) → gffread `-w` → `RENAME_FASTA_HEADERS` writes the
renamed FASTA **and** `<prefix>.id_map.tsv` (renamed header → original gffread seqname,
which is the GFF3 `ID=` verbatim) → TD2 keeps the non-coding records (lncRNA only) → that
final FASTA is emitted as `spliced_fasta` **and** handed, with the candidates GFF3 and the id
map, to `FILTER_ANNOTATION` (`filter_gff_by_id.py` keep mode).

Keep mode matches by **exact id** — never the `normalise()` prefix stripping of drop mode,
which disagreed with `RENAME_FASTA_HEADERS`' colon split for ids like `transcript:foo:bar`
and could leave a TD2-dropped transcript in the GFF3. It exits 1 if a renamed header maps
twice, a FASTA record is missing from the map, or a mapped id is not a transcript of the
GFF3 — so on success the GFF3's transcripts **are** the FASTA's records. The joins feeding
it use `failOnMismatch: true` for the same reason.

Consequences: the reference is longest-isoform **and TD2-filtered**, deliberately — the
projected models are TD2-filtered too, so an unfiltered reference counted coding lncRNA
the pipeline can never report as false negatives. Transcript-level Sn is now ≈ per gene,
precision can dip when a source's longest isoform matches a *non-longest* target isoform,
and numbers are not comparable to pre-plans/20 runs. `TD2_NONCODING.out.coding_ids` is
still emitted but no longer consumed. With `--longest_isoform false` the projected FASTA
and the filtered GFF3 are byte-identical to the pre-plans/20 pipeline (verified); only the
reference (now TD2-filtered) differs.

### The merge/split contract

Merging every source into one alignment is lossless because the species is already carried
end to end, and **nothing in the chain may break that chain of custody**:

`modules/local/rename_fasta_headers.nf` writes `>tid|type|species` → minimap2 puts that in
the BAM QNAME → `bin/bam_to_gff.sh` stamps it on every `gene` and `transcript` row as
`source=<species>` → `bin/gff_by_source.py` keys on that attribute (child `exon`/CDS rows
inherit it through `Parent`), writing one file per source (`--name-template` with a
`{source}` placeholder; its subset mode went with the top-N consensus, plans/23). The same
header species field is what MINIMAP2_ALIGN's `exclude_self` filter reads, so self-exclusion
and the split can never disagree on who a record belongs to.

Four things that will bite:

- **The split filename template and the re-keying arithmetic must agree.** `ext.args` for
  `SPLIT_GFF_BY_SOURCE` builds `<target>.from_{source}.<gtype>.projected.gff3`. Since
  plans/18 that per-source list is never re-keyed in `projection.nf` itself (it's folded,
  unexploded, into `GFF_STATS_PROJECTED_BATCH`'s input list) — the re-keying now happens one
  step further downstream, in `benchmarking.nf`, which strips the analogous
  `<target>.from_` / `.<gtype>[.decoy].gffcompare.stats` prefix/suffix *by length* off
  `GFFCOMPARE_BATCH`'s output filenames (not by regex — species names are full of `_` and
  digits) to recover each source's identity for the per-source accuracy. Change the split
  template and that re-key breaks too, with an out-of-range substring — it does not silently
  mismatch, which is the point.
- **Transcript IDs must be unique across the merged file.** Per-source GFF3s were separate
  files, so a shared id never collided; in one file a duplicate `ID=` makes `gffread -w` and
  gffcompare (benchmarking, curation) silently misbehave. Source transcript ids are **not** guaranteed
  unique across species — seen for real, two *Drosophila* species both used the generic id
  `lnc_RNA1412` — so `bam_to_gff.sh` builds the projected `ID=`/`Parent=`/`gene-` prefix as
  `<tid>|<source>`, not bare `<tid>`, making it globally unique by construction rather than
  assuming upstream ids never collide. `gff_by_source.py` still asserts uniqueness and exits
  1 naming the id and both sources if that guarantee is ever broken — it should never fire
  now, but fails loudly rather than corrupting the per-source split if it does. One knock-on: TD2's
  coding-potential drop-list (`bin/td2_coding_filter.py`) must NOT strip this qualified id
  back to a bare tid for the *projected* stage the way it correctly does for the *input*
  stage (single source per file there, so bare tid is unambiguous) — doing so would let one
  source's coding call drop an unrelated same-named transcript from a different source.
- **One source yields a bare Path, not a List.** The `path("*.gff3")` glob of
  `SPLIT_GFF_BY_SOURCE` is a List for ≥ 2 sources and a bare Path for one (common now that a
  `both` target loses its self file); `projection.nf` normalises with
  `gffs instanceof List ? gffs : [gffs]`.
- **A target with nothing projected onto it is real, not broken — and `SPLIT_GFF_BY_SOURCE`'s
  output is `optional: true` because of it.** If `allModels.raw.gff3` has zero records
  (nothing from any source survived alignment + TD2 filtering onto this target — seen for
  real on divergent targets in a large all-vs-all run, and — since plans/23 — a `both`
  target whose only source is itself), `gff_by_source.py` writes zero files and exits 0 by
  design. `projection.nf` handles it with `join(..., remainder: true)` rather than a plain inner join, which would otherwise make the
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
    `'.*:GFF_STATS(_TARGET|_PROJECTED(_BATCH)?)?$'` — keep it in sync if a new
    `GFF_STATS` alias is added. **`GFF_STATS_PROJECTED_BATCH`** (plans/18) is the one
    exception to the `ubuntu:22.04` base: it also needs python3 (to run
    `gff_stats_to_mqc.py` in the same task, see Reporting conventions), so it runs on
    `python:3.11` instead. **Not yet verified** to meet the glibc ≥ 2.34 floor — check
    with `docker run --rm python:3.11 ldd --version` before relying on it.
  - It cannot compute genome coverage (no genome-size input).
  - It resolves a transcript's gene via `Parent`/`gene`/`Gene` only — so a gffread-derived
    GFF3, which carries the gene in `geneID=` and emits no gene features, reports 0 genes
    (none is produced since plans/23; the curated GFF3s carry real gene rows).
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
`modules/local/{samtools_to_mqc,gffcompare_to_scatter}.nf`, which reuse the MultiQC
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
   `*_transcript_filter_mqc.tsv`,
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
