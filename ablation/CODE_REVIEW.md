# FOMO codebase review — a maintainability "ablation"

*2026-10-15 · branch `ablation` · code at `2bf111c` · scope: everything tracked except `plans/`,
test data and `ablation/`*

The question asked: what in this codebase is over-complicated, duplicated, inconsistently named,
needlessly nf-core-shaped, run-specific or over-verbose — and what strategy keeps it
maintainable and easy to extend? The companion [`ABLATION_REPORT.md`](ABLATION_REPORT.md) asked
what each *scientific* component is worth; this one asks what each *line of code* is worth.

**Method.** Four independent read-only review passes (Nextflow workflows/config, modules, Python,
shell/assets/docs), every claim checked against the code, plus:

- automated metrics — comment density and "debt markers" (`scripts/code_metrics.py`), ruff, `nextflow lint`;
- git hotspots;
- three experiments: a stub run, the strict v2 parser, and the batch-loop exit codes.

The high-impact claims were re-verified by hand. The verification log is at the end.

---

## Summary

The pipeline is carefully *reasoned*. Almost every odd-looking line has a recorded reason, and
the data contracts ("fail loudly", one filtered annotation, merge/split custody) are sound.
The maintainability problems come from **where that reasoning is stored**. Three root causes
explain most findings:

1. **Identity is encoded in strings and recovered by parsing, instead of carried in the meta
   map.** The track name is assembled inline 35 times. The same filename is parsed by three
   different parsers (Groovy, bash, Python) with different rules. MultiQC *display* labels
   (`1_raw`, `2_lncRNA`) are parsed back as *data keys* by the run summary. Adding a feature
   type touches 8+ files in three languages, and a mismatch usually drops rows silently rather
   than failing.
2. **Prose stands in for structure and tests.**
   - The Nextflow layer has almost one comment line per code line (0.85).
   - `CLAUDE.md` is 514 lines and reads as a changelog.
   - ≥ 20 comments now contradict the code, including user-facing `--help` text.
   - Comments cite 23+ specific runs, job IDs and clades.
   - 13 of 14 scripts and both core awk programs have no tests; there is no CI.
   - The stub run, which every module carries a stub for, does not complete.

   Invariants that a test could enforce are instead re-explained in up to six places, which
   have already drifted apart.
3. **nf-core form without nf-core tooling.**
   - Two "nf-core" modules are local forks; `nf-core modules update` on `minimap2/align`
     would silently bring self-projection back.
   - Two nf-core modules are never used.
   - 25 version outputs are emitted and never collected.
   - Labels exist only to be overridden, and the conda profile cannot run 18 of 21 local
     modules.

   The boilerplate costs maintenance without delivering the update and lint benefits it
   exists for.

**Strategy in one line:** move identity into meta and a registry, move invariants into tests, keep
`modules/nf-core` for unmodified tool-managed modules only, and cut prose to present-tense "why".
Do it in five phases, **safety net first** (§5).

| what | today | target |
|:---|---:|---:|
| comment : code, Nextflow workflows + config | 0.83 / 0.88 | ≤ 0.4 |
| `CLAUDE.md` | 514 lines | ~150–200 |
| copies of the track-name expression | 35 in 13 files | 1 (`meta.track`) |
| parsers of the per-source filename grammar | 3 (+1 partial) | 1, or none (manifest) |
| scripts / awk programs with tests | 1 / 0 of 14 / 2 | all pure helpers + fixtures |
| `nextflow lint` errors (strict parser) | 6 (pipeline fails to compile under v2) | 0 |
| Python functions over ruff C901 (complexity 10) | 10 (max 38) | 0 above 15 |

---

## 1. Metrics

### Comment density and documentation debt

*`scripts/code_metrics.py`; comment lines per code line*

| layer | code | comments | ratio | `plans/NN` refs | history phrases | run / incident refs |
|:---|---:|---:|---:|---:|---:|---:|
| Nextflow config (`nextflow.config`, `conf/`) | 415 | 366 | **0.88** | 27 | 13 | 25 |
| Nextflow workflows + subworkflows | 452 | 373 | **0.83** | 21 | 2 | 2 |
| `modules/local` | 728 | 202 | 0.28 | 12 | 13 | 4 |
| `modules/nf-core` | 461 | 41 | 0.09 | 2 | 0 | 0 |
| Python `bin/` | 2,467 | 492 | 0.20 | 23 | 9 | 4 |
| shell `bin/` | 341 | 102 | 0.30 | 2 | 1 | 5 |
| `CLAUDE.md` | 514 lines | — | — | 41 | — | — |

The phrase and reference counts come from a regex and are lower bounds; the reviews found more
by reading (e.g. ≥ 41 history-narrating lines in the Nextflow layer alone).

The worst single files are:

- `conf/crg.config` 1.81
- `conf/base.config` 1.53
- `modules/local/gff_stats.nf` 1.50
- `projection.nf` 1.10
- `benchmarking.nf` 1.10
- `td2_predict.nf` 1.13

The comment density is inverted relative to the complexity. The densest code
(`curate_models.py`, McCabe 37–38) has 0.13. The simplest declarative files (resource configs)
have 1.5–1.8.

### Static analysis

- **ruff** (E, F, B, C90, N, UP, SIM, PL): 283 findings.
  - No unused imports or names.
  - 10 functions over complexity 10:
    - `curate_models.main` 37 (163 statements);
    - `curate_models.read_reference` 38;
    - `filter_gff_by_id.main` 21;
    - `curate_models.read_models` 20;
    - `write_merged` 19;
    - `gff_by_source.main` 17.
  - The test's `main` has 187 statements.
- **`nextflow lint`** (Nextflow 25.10):
  - 6 errors: `reporting.nf:47` and `run_summary.nf:34` fail to parse
    (`{ _meta, p -> p instanceof List ? p : [p] }`), which cascades into 4 more.
  - 15 warnings: deprecated `Channel`, implicit `it`, unused closure parameters.
  - Confirmed by running with `NXF_SYNTAX_PARSER=v2`: *"Script compilation failed"*. The
    pipeline will stop compiling when the strict parser becomes the default.

### Change hotspots

Of 72 commits, the most-changed files are:

| file | commits |
|:---|---:|
| `conf/modules.config` | 28 |
| `nextflow.config` | 27 |
| `CLAUDE.md` | 24 |
| `workflows/fomo.nf` | 16 |
| `projection.nf` | 14 |
| `preprocessing.nf` | 12 |
| `benchmarking.nf` | 12 |

These are exactly the files with the highest comment density and the most duplicated
string grammar. They change the most and are the costliest to change.

---

## 2. Code ablation: what can be deleted with no change in behaviour

Everything below is dead, unused or redundant today. Removing it changes no output (the two
starred rows change only duplicate report rows).

| item | evidence | size |
|:---|:---|---:|
| `modules/nf-core/gunzip`, `modules/nf-core/gffcompare` | included nowhere. They are referenced only implicitly: `gffcompare_batch.nf:18` borrows the conda env path, and `merge_fasta.nf` borrows the gunzip image just to run `cat` | 125 lines |
| `modules/nf-core/.DS_Store` | tracked despite `.gitignore` | 1 file |
| unused `minimap2/align` branches | FOMO always passes `true, 'csi', false, false`. The PAF, BAM-input and `samtools reset` paths, `args3`/`args4` and the else-branch never run | ~35 lines |
| unused `gffread` outputs and branches | all 3 call sites pass `-w` with a FASTA. The gtf/gff3/bed outputs, extension inference and the no-FASTA branch are dead | ~25 lines |
| per-source re-key in `benchmarking.nf:77-86` | the `id` it reconstructs is read by nothing: REPORTING uses `target_id`, RUN_SUMMARY drops the meta, and the scatter re-parses filenames in Python. Replace it with `.transpose()` and delete the 3 coupling comments | ~20 lines |
| string-key join `projection.nf:198-206` | both sides carry the same meta; `join(..., remainder: true)` on meta works, and `td2_noncoding.nf` already does this | ~15 lines |
| `GFF_STATS_TARGET` + `_TO_MQC` for `both` species* | recomputes stats of the same two files PREPROCESSING already processed. Only "raw stats of a pure target with gff3" is unique | 2 aliases, ~40 config lines |
| unused config | `withLabel: 'big_mem' / 'small'` (no module uses them); a 19-line commented GPU block with 3 TODOs; MINIMAP2_INDEX `ext.prefix` (never read); the dead selector alternatives `GFFCOMPARE\|` and `_PROJECTED`; six `ext.prefix` that restate module defaults | ~60 lines |
| `versions.yml` filters in 8 `publishDir` blocks | no published process writes `versions.yml` | 8 closures |
| emits nobody consumes | `PREPROCESSING.intergenic_bed/decoy_gff3`, `PROJECTION.bam/index`, `CURATION.curated/merged/report`, the REPORTING/RUN_SUMMARY emits, each with a cardinality comment to keep accurate | ~15 lines |
| dead Python | `relocate_loci.build_shape`, `Feature.attrs_raw`, `fomo_stats.Identity.is_aggregate/is_self` (re-implemented inline), `read_models`' discarded second return, the `top\d+` filename alternative, `gff_stats_to_mqc --section-id` (only logged) | ~40 lines |
| `empty_ids.txt` workaround | make the drop list an optional `[]` input | 1 asset + 6 lines |
| `assets/samplesheet_full.csv` | used nowhere; its paths do not exist (`/Users/fzanarello/...`, `assets/annotation_downloads/`) | 1 file |
| test species *L. phlaeas*, *L. virgaureae* | in no samplesheet | 26 MB |
| `--include_single_exon`* | no downstream effect: benchmarking runs `-M` and curation drops single-exon models (ablation B5) | 1 param |
| stale, history and incident comments | §3.5 | **~450 lines** if the Nextflow layer is brought to 0.3 |
| `CLAUDE.md` | §3.5 | **~330 lines** |

That is roughly **1,000+ lines (≈ 15 % of the non-test source plus docs)** removable before any
real refactoring.

---

## 3. Findings

Severity: **H** = blocks safe extension or hides errors · **M** = recurring cost · **L** = polish.

### 3.1 Over-complicated paths and functions

- **P1 · H · identity in filenames.** Model identity `<target>.from_<source>.<gtype>` is written
  into filenames and parsed back three ways:
  - a Groovy length-slice, `benchmarking.nf:80-82`;
  - bash `${id%%.*}`, `gff_stats_projected_batch.nf:76-77`, whose comment says ids "never
    contain a literal '.'";
  - Python `FNAME_RE`, `fomo_stats.py:41-43`, hard-coding `lnc_RNA|mRNA`;
  - plus a partial `split(".from_")` in `run_summary_tables.py:265`.

  `assets/schema_input.json:12` allows any `^\S+$` species, so `.` truncates the bash parser,
  and `|` or `'` would break the `>tid|type|species` header and `awk -v s='…'`.

  *Fix:* identity lives in meta. Batch tasks write a small manifest (file, source, target,
  track) instead of being parsed. Tighten the species pattern to `^[A-Za-z0-9_-]+$`.
- **P2 · H · display labels used as data.** The MultiQC row-ordering prefixes are built in
  Groovy in five ladders with two numbering schemes:
  - `1_raw/2_lncRNA/3_decoy_lncRNA/4_mRNA/5_decoy_mRNA` at `modules.config:87,139`;
  - `1_lncRNA/2_decoy_lncRNA/…` at `modules.config:268` and `gff_stats_projected_batch.nf:54`.

  They are parsed back as keys by `run_summary_tables.py` (`PROJECTION_TOKEN`:49, hard-coded
  `"1_raw"`/`"2_lncRNA"` at :374-375, :428-429, :549, :552, `"1_lncRNA"` at :773). A
  mismatch drops summary rows silently.

  *Fix:* keep `order` and `label` in meta, set from one registry. The summary reads
  machine-readable stats rather than presentation rows (P4).
- **P3 · H · MultiQC suffix routing is a three-sided contract with no check.** It spans:
  - the producers' `*_xxx_mqc.tsv` suffixes (7 modules and 2 scripts);
  - `sp` patterns in `assets/multiqc/sections.yml:153-173`;
  - `run_summary_tables.py:144-155` `pick("_xxx_mqc.tsv")`;
  - plus section ids against `report_section_order` keys.

  A typo makes `pick()` return `[]` and leaves a section without an order — silently.

  *Fix:* one section registry (YAML or a Python module) read by the summary, plus a test that
  every `sp` suffix has a producer and an order entry.
- **P4 · M · the run summary scrapes presentation files.** `RUN_SUMMARY_TABLES` reconstructs
  every number from MultiQC custom-content TSVs:
  - it rebuilds the same indexes repeatedly: `accuracy()` 3×, `annotation_index()` 3×,
    `projection_index`/`td2_index` 2× each;
  - each `section_*` builder (32–102 lines) mixes computing, warnings, TSV, JSON and an
    out-parameter dict.

  This is the task recorded hitting its time limit on the 531×1,260 run (`crg.config:69-80`).

  *Fix:* `functools.cached_property` on `Inputs`, and pure `compute_<section>()` plus
  `emit_<section>()`. Longer term, have producers emit data JSON, render MultiQC from it, and
  aggregate from the data.
- **P5 · M · the config include order creates the "precedence trap".** `nextflow.config`
  includes `base.config` *after* `profiles{}` (lines 75–123), so base selectors override
  profile selectors. About 20 comment lines in `base.config`/`crg.config` exist only to warn
  about it.

  *Fix:* use the nf-core order (base → profiles → modules) and delete the warnings. A reviewer
  confirmed with `nextflow config` that this restores profile precedence.
- **P6 · M · fully qualified selectors make the call tree rigid.**
  - 34 of 45 `withName` selectors are `FOMO:SUB:PROC`.
  - This is the only reason for CLAUDE.md's "never nest subworkflows" rule, for `curation.nf:9-11`,
    and for making RUN_SUMMARY a separate top-level subworkflow.
  - Anchoring is inconsistent (15 anchored, 27 unanchored qualified selectors), although
    `modules.config:61-63` warns that unanchored ones leak.

  *Fix:* use `'.*:NAME$'` everywhere. Disambiguate the two TD2 instances by alias
  (`TD2_CODING_FILTER_INPUT/_PROJECTED`) or by a meta key.
- **P7 · M · nine copies of the same channel shape.** The pattern is re-key by id →
  `combine(by:0)` → multiMap/strip the key:
  - `preprocessing.nf:56,93,105,120,132,203`
  - `projection.nf:140`
  - `benchmarking.nf:53`
  - `curation.nf:30`

  `FASTA_SIZES.out.sizes.map{…}` is built twice (`preprocessing.nf:96,108`).

  *Fix:* keep id-keyed channels through PREPROCESSING, or use one `pairById()` helper in a
  `lib/`/functions file.
- **P8 · H · god functions in Python.**
  - `curate_models.main` has 228 lines and does nine jobs (CLI, track inference, filtering,
    reference tagging, chain support, discard report, gene building and hashing, novel
    stats, four writers).
  - `read_reference` has 137 lines with 13 mutable locals shared by 3 nested closures, GTF
    and GFF3 interleaved.
  - The test reaches one branch by monkey-patching `gene_hash`, `sys.argv` and `sys.stderr`.

  *Fix:* split into `resolve_track / filter_models / tag_ref_hits / chain_support /
  build_genes / novel_report / write_*`, and use a `_RefBuilder` class with
  `add_gtf_row`/`add_gff3_row`.
- **P9 · M · core selection logic in a two-pass awk.** `bin/filter_transcript.sh`, the
  pipeline's *only* isoform selection point, is ~125 lines of awk that reads the file twice.
  Its two attribute parsers are stricter than the Python ones (no whitespace trimming). The
  double read forces the scratch-decompression boilerplate in `filter_transcript.nf:30-39`.

  *Fix:* rewrite it in Python on the shared GFF module (reads `.gz` natively, no scratch dir),
  with fixtures.
- **P10 · L · three gating mechanisms for optional features.**
  - decoys use `if (params.include_decoy)` in a subworkflow;
  - curation uses `ext.when`;
  - mRNA uses the feature-type list.

  *Fix:* pick one pattern and document it once.

### 3.2 Duplicated functions and logic

- **D1 · H · GFF parsing: six Python variants plus awk, with different semantics.**
  - `attr()` is byte-identical in `filter_gff_by_id.py:49` and `gff_by_source.py:48`.
  - `TX_TYPES` is copied verbatim in the same two files ("Mirrors…").
  - `curate_models` defines "transcript" the other way round (`GENE_LEVEL` + `NOT_TRANSCRIPT`,
    case-sensitive), so an unlisted type is a transcript to one script and not to another.
  - `relocate_loci.parse_gff3_attrs` keeps the *first* duplicate key; `curate_models.gff3_attrs`
    keeps the *last*.
  - The "skip #, split tab, need 9 columns" loop appears ~10 times.

  *Fix:* `bin/fomo_gff.py` (`open_text`, `iter_gff`, `parse_attrs`, `FEATURE_VOCAB`),
  imported the way `fomo_stats.py` already is.
- **D2 · M · MultiQC helpers.**
  - 5 TSV writers.
  - 2 `fmt()` functions that disagree on NaN (`"nan"` vs `"NA"`).
  - The custom-content JSON skeleton repeated 5×.
  - Percent computed 3 ways.

  *Fix:* `bin/fomo_mqc.py`.
- **D3 · M · the two batch modules are one scaffold copied.** `gffcompare_batch.nf` and
  `gff_stats_projected_batch.nf` share:
  - tag and label;
  - the aggregate-symlink trick (4 copies across script and stub);
  - the `run_one` / `export -f` / `xargs -P` loop;
  - 30+ lines each narrating the "old ext.prefix".

  *Fix:* one `bin/for_each_projected.sh` with strict mode (see R1).
- **D4 · M · GFF statistics is three modules, two images and four near-identical config blocks.**
  - `GFF_STATS` + `GFF_STATS_TO_MQC` (2 tasks per file) against `GFF_STATS_PROJECTED_BATCH`
    (both in one task).
  - The same vendored binary runs on `ubuntu:22.04` in one and `python:3.11` in the other.
  - The kind→suffix closure is copied 4× (`modules.config:65,78,115,127`).
  - `FILTER_ANNOTATION` and `FILTER_GFF_BY_ID` are two modules around one script.
  - Scratch staging is copied (`gffread`, `filter_transcript`).
  - The Groovy `.gz` branching appears in 4 modules, where `gzip -cdf` would do.
- **D5 · M · the same explanation in up to six places.** These have already drifted:
  - **Trans-spliced genes:** listed as `nad1/nad2/rps12` in CLAUDE.md, `nad1/nad2/nad5` in
    `filter_transcript.sh`, and `nad1/nad2` in `sections.yml`.
  - Unique IDs and the *Drosophila* case: CLAUDE.md, `bam_to_gff.sh`, `gff_by_source.py`.
  - The gzip/faidx warning: CLAUDE.md, `preprocessing.nf`, `projection.nf`, `fasta_sizes.nf`.
  - "combine spreads a List": CLAUDE.md, `fomo.nf`, `preprocessing.nf`, `reporting.nf`.
  - The `both` role: README, CLAUDE.md, both schemas, `fomo.nf`, `benchmarking.nf`,
    `test.config`.
  - The output tree: CLAUDE.md, README, `modules.config`.
- **D6 · L · repeated literals.**
  - TD2 args 3× (`td2_predict.nf:38`, `modules.config:39,233`).
  - The default 0.5 three times and 1000 twice (script defaults vs `nextflow.config`).
  - `LOCATIONS` in two scripts in different orders.
  - Curation metric names as string literals in `curate_models`, `curation_to_mqc` and
    `run_summary_tables`, read with `.get(k, 0)`, so a renamed metric silently becomes 0.

### 3.3 Inconsistent naming

- **N1 · H · one concept, four spellings.** The (feature type, decoy) track is spelled:
  - `lnc_RNA.decoy` (filenames, `--track`);
  - `lnc_RNA_decoy` (`meta.gtype`, set at `projection.nf:50` and used *once*);
  - `decoy_lncRNA` (MultiQC tokens, curated ID tags);
  - `lnc` (branch names, `EXTRACT_PROJECTED_LNC`).

  The inline `${meta.feature_type}${meta.decoy ? '.decoy' : ''}` appears **35×** in 13 files,
  including four exact copies of `"${meta.target_id}.allModels.<track>.raw"` (`modules.config:215,246,261,265`).

  *Fix:* stamp `meta.track` (dot form) and `meta.label` once, in PREPROCESSING. A meta key is
  the reliable option because config closures cannot call `lib/` helpers dependably.
- **N2 · M · `meta.id` is overloaded.** It goes source species → track → the literal
  `'allModels'` → `'curated'` → target.
  - Consequence: the task tags for MINIMAP2_ALIGN, SAMTOOLS_STATS, FILTER_ALLMODELS and the
    projected TD2 are identical for every target (`allModels`), so traces and logs cannot tell
    targets apart (seen in the Part B traces of the ablation study).
  - *Fix:* explicit `source` / `target` / `track` keys, `id` never repurposed, and one tag rule.
- **N3 · M · emit and alias names.**
  - MultiQC emits are called `tsv` (3 modules), `mqc` (5) and `json`.
  - `TD2_NONCODING` emits `mqc` where every other subworkflow emits `mqc_files`.
  - Alias styles mix verbs and objects: `FILTER_GFF_BY_ID→FILTER_ALLMODELS`,
    `GFF_BY_SOURCE→SPLIT_GFF_BY_SOURCE`, `MERGE_FASTA→MERGE_SOURCE_FASTA`, `TD2_NONCODING_PROJ`.
  - `allModels` vs `allmodels_raw` vs `allModels_raw` vs `allsources`.
  - One include keeps `main.nf` in its path.
- **N4 · M · Python CLI conventions.**
  - Input flags: `--gff`, `--input-gff`, `--in-fasta`, `--input`, `--json`, a positional, `--mqc-dir`.
  - Sample flags: `--name` vs `--sample`.
  - Feature flags: `--feature-type` vs `--track` vs `--feature lncRNA`.
  - Booleans: `--include-mrna true` vs `store_true`.
  - Four logging dialects, and exit codes 1/2 or `None`.
  - `curate_models.py` uses single quotes, no type hints and 48 % of locals ≤ 2 characters
    (`R`, `obs`, `dis`, `nov`; `m` both a regex match and a Model). The other scripts use
    double quotes, typing, and the old `typing.Dict/List` (105 UP006).
- **N5 · L · small inconsistencies.**
  - `gff-feature-stats` (upstream kebab-case; say so).
  - The `*.neflow_command*.txt` typo is the reason for three `.gitignore` patterns.
  - The `relocate_loci.py` help says decoys go into the "Target intergenic intervals"; they go
    into the *source's*.

### 3.4 Unnecessary nf-core complexity

- **C1 · H · forked modules under `modules/nf-core`.**
  - `minimap2/align` holds core FOMO behaviour: `exclude_self` (no self-projection) and a
    no-pipe rewrite.
  - `gffread` holds the scratch staging.
  - There is no `modules.json`, no patch file, no `meta.yml` and no tests. **Running
    `nf-core modules update` would silently restore self-projection.**
  - `exclude_self` lives only in the BAM branch, so the module called with
    `bam_format=false` would self-project without any error.
  - Upstream `gffread` now handles `.gz` itself, so the fork no longer tracks upstream.

  *Fix:* move both to `modules/local` and cut them to what FOMO calls:
  - `minimap2_align_splice`: `tuple(meta, fasta, mmi)` → exclude-self, minimap2, sort with CSI;
  - `extract_transcripts`: stage → `gffread -w`.
- **C2 · M · keep nf-core modules only where the tooling is used.** Four modules are unchanged
  apart from container or version: `minimap2/index`, `multiqc`, `samtools/stats`,
  `seqkit/stats`. `bedtools/complement` also qualifies. Install these with `nf-core modules
  install` (which creates `modules.json`) so `update` and `lint` work and drift becomes
  visible. Delete `gunzip` and `gffcompare` (§2).
- **C3 · M · version reporting is half-wired.**
  - 25 `topic: versions` outputs; **no `channel.topic('versions')` consumer exists**.
  - `FILTER_TRANSCRIPT` alone still writes a legacy `versions.yml`.
  - 6 local modules report nothing.
  - 4 `awk --version` evals return an error string with Ubuntu's mawk.
  - The manifest floor (`!>=24.10.0`) is probably below the release where topic channels
    stopped needing a preview flag.

  *Fix:* either wire one topic collector into MultiQC (~15 lines) and fill the gaps, or delete
  all of it. Do not keep paying for unread output.
- **C4 · M · ext knobs, constant inputs and placeholder labels.**
  - Behaviour is set through `ext.role/section_id/genes_output/stage/sample_name`.
  - `RELOCATE_LOCI` alone reads `params.*` inside the module.
  - `when:` appears in 30 modules but `ext.when` is used once.
  - Constant `val` inputs and unused `meta2`/`ref_meta` inputs.
  - Labels that are immediately overridden: `process_high` → 4 cpu / 4 GB for both of its
    users; `process_medium`'s 4 cpu → 2 for its only user.
  - MINIMAP2_INDEX is labelled `process_low` (4 GB) while comments say a real genome needs
    25–30 GB.
- **C5 · M · containers.**
  - Five image families for code that is mostly stdlib Python and coreutils.
  - Floating tags: `python:3.11` in 8 modules, `ubuntu:22.04` in 6.
  - The MultiQC image is used to run stdlib Python scripts, which ties them to MultiQC bumps.
  - Two Wave URL styles.
  - `-profile conda` cannot run 18 of 21 local modules (no `conda` directive).

  *Fix:* three images pinned by digest (python-stdlib, coreutils/awk, curate =
  python+gffcompare+gffread), set by label in config. Then drop the conda profile, or give every
  module an `environment.yml`.

**Verdict on the boilerplate:**

- For unmodified, tool-managed nf-core modules, keep the full pattern (meta, `ext.args`,
  `ext.prefix`, topic versions).
- For local one-liners, keep `meta` for routing and publishing. Keep `ext.args` only where
  config really sets it. Replace `ext.prefix` with defaults built from `meta.track`, and make
  the versions decision once.

### 3.5 Comments: run-specific, historical, stale, verbose

- **V1 · H · stale comments that contradict the code (≥ 20), including user-facing text.**

  User-facing:
  - **`nextflow_schema.json:103`** (`--help`) says curation reports "one representative per
    gene (most splice junctions, then most species…)". It reports one transcript per chain,
    t1 = most species, and the text ends with "See plans/21".
  - `:22` mentions the removed "top-N consensus".
  - `:120` says TD2 is "per-source"; it is pooled per target.
  - `sections.yml:65,94` describe per-source alignment rows and "consensus".
  - `main.yml:10` hard-codes "at least two species" instead of `--curate_min_species`.
  - **`nextflow.config:126-127`** promises timestamped `pipeline_info` reports "so every run
    keeps its own history". All four use fixed names with `overwrite = true`, **so the traces
    used for resource right-sizing are overwritten by every run.**

  Internal:
  - `fomo.nf:83-84` says "the schema rejects a gff3 on a pure 'target'" (false since plans/22).
  - `modules.config:110` says the target reference is deliberately *not* TD2-filtered;
    `benchmarking.nf:29` says it is, on purpose.
  - Three comments place the re-key in `projection.nf`; it is in `benchmarking.nf`.
  - Names of removed things are still cited: `ANNOT_ORDER`, `FILTER_TARGET`,
    `GFF_STATS_PROJECTED`, `GFFCOMPARE_COMBINE`, "top-3 consensus".
  - `td2_predict.nf:70-74` calls zero ORFs "the norm" (self-projection is gone).
  - `test.config` says "4 target×source pairs" and "48 alignments" (actually 3 and 8).
  - `crg.config:18` cites 12 GB / 14 GB caps that no longer exist.
  - Drifted line references: `crg.config:43` "base.config:80-92", `crg.config:73`
    "run_summary_tables.py:138".
  - In Python: `run_summary_tables.py:215`'s docstring states the *opposite* merge precedence
    of the code (R2). `run_gffcompare`/`exact_chains` docstrings give the wrong return types.
- **V2 · M · run and incident references in code (23+).**
  - `crg.config:50` "OOM'd for real on the mollusca run (jobId 28546390: ReqMem=24G,
    MaxRSS=25,163,300K…)"
  - `crg.config:71-77` "Lepidoptera run (S=531, T=1260) … job 28348885"
  - `crg.config:123` "run `agitated_legentil`"
  - `base.config:29-31` "was 12.GB — worst observed peak 0.23 GB (GFFREAD)"
  - `modules.config:192` *Vampyroteuthis infernalis*
  - `projection.nf:85` "a tick chromosome at 560,970,318 bp"
  - `curate_models.nf:3-4` "82 targets (191 k models from 176 species…)"
  - `td2_predict.nf:48` "CRG has `mount home = yes`"
  - `fomo.nf:28` "Verified the hard way."

  The rule each one justifies is worth one line; the anecdote belongs in the commit message or
  the plan.
- **V3 · M · history narration.**
  - ≥ 41 lines in the Nextflow layer: "used to", "no longer", "the old ext.prefix", "un-batched",
    "since plans/NN".
  - `gff_stats_projected_batch.nf:36-53` is an 18-line story of a prefix that no longer exists.
  - `nextflow.config:45-73` has 29 lines about where publishDir *used to* live and why conda
    *used to* be global.

  Narration describes a change, not the code, so every refactor makes it wrong — which is
  why V1 exists.
- **V4 · M · verbosity.**

  | comment lines | code lines | where |
  |---:|---:|:---|
  | 12 | 1 | `fomo.nf:17-28` |
  | 25 | 9 | `projection.nf:173-197` |
  | 7 | 1 argument (`'csi'`) | `projection.nf:83-89` |
  | 20 | 2 settings | `base.config:59-78` |
  | 16 | 2 settings | `crg.config:48-63` |

  - `gff_stats.nf:5-34` restates CLAUDE.md "Key Tools" almost word for word.
  - Some comments describe things that are *absent*: `modules.config:279-283` (a block for
    `GFF_STATS_PROJECTED_BATCH` that is not there), `:286-288` (a `FILTER_TARGET` that no
    longer exists).
- **V5 · M · `plans/NN` as justification.**
  - 92 references in code, 41 in CLAUDE.md, and one in user-facing `--help`.
  - Only 9 of 30 plans carry a status. `plans/09` is superseded but unmarked. `plans/25` says
    "not yet committed" but is committed.

  A reference such as "(plans/18)" sends the reader to a document that mixes the decision,
  rejected alternatives and run results.
- **V6 · M · `CLAUDE.md` is a changelog.** It is 514 lines and is organised by when features
  were added:
  - the gunzip/I-O section sits inside "Feature tracks";
  - the legacy AGAT table contradicts the "AGAT is no longer used" note further down;
  - it carries incident numbers ("~600 GB per run on the 15-cephalopod set") and an open TODO
    (the glibc floor is "Not yet verified").

  It is the single biggest duplication source (D5).
- **V7 · L · the README is not a user entry point.**
  - `nextflow run . -profile test` enables no container engine (only `crg` does).
  - It has no requirements, samplesheet example, parameter table or contributor notes.
  - It sends readers to CLAUDE.md for architecture.

### 3.6 Robustness and standard practice

- **R1 · H · the batch modules mask failures and defeat retry-on-OOM.** In the `xargs … bash -c
  'run_one'` worker:
  - a failing non-final command is ignored, because the child shell has no `-e`;
  - a killed worker makes `xargs` exit **123**, not 137, so `base.config`'s retry rule
    (`130..145, 104`) never fires.

  Verified: `false; echo` gives exit 0, and `kill -9` gives 123. This is the exact masking that
  the minimap2 rewrite fixed, still present in `GFFCOMPARE_BATCH` and
  `GFF_STATS_PROJECTED_BATCH`.

  *Fix:* `set -euo pipefail` in the worker, plus 123/125 in those processes' retry list (or
  `parallel --halt now,fail=1`).
- **R2 · M · latent bug in `run_summary_tables.curation_index` (:214-222).** It merges the
  *observed* and *novel* rows with `dict.update`. Both have a `genes` column with different
  meanings (pre-filter vs curated). The result is correct only because the novel rows are
  concatenated last. The docstring says the opposite precedence, so "fixing" the code to match
  the docstring would report pre-filter genes as curated.

  *Fix:* namespace the columns, or keep two dicts.
- **R3 · M · the chain of custody is not enforced at its entry point.** `bam_to_gff.sh:77-80`
  silently sets the source to `.` when a read name has fewer than three `|` fields.
  `gff_by_source.py` then writes `<target>.from_..<track>…` and exits 0.

  *Fix:* `np != 3` → error.
- **R4 · M · strict-syntax and deprecation debt.**
  - The pipeline does not compile under the v2 parser (§1).
  - It uses deprecated `Channel.` (5×) and implicit `it` (3×).
  - The config has tabs in 3 files.
- **R5 · M · error handling and logging.**
  - `sys.exit` inside library functions: `curate_models.die`, 7× in `run_summary.Inputs`.
  - Silent fallbacks in the same class.
  - Exceptions swallowed in `samtools_to_mqc:74`.
  - A malformed PSAURON score is silently treated as non-coding (`td2_coding_filter:86`).
  - `assert` used for a runtime check (`curate_models:856`).
  - A temp dir that leaks on error (`curate_models:844`).

  *Fix:* `FomoError` raised, and `main()` logs and exits 1.
- **R6 · L · `relocate_loci.py`.**
  - It uses the global `random` state, so determinism depends on call order. Use
    `rng = random.Random(seed)`.
  - It recomputes chromosome weights per model, which is O(models × intervals); with
    `--decoy_cap 0` that becomes 10⁴–10⁵ squared.
- **R7 · L · site specifics in the shared repo.**
  - `crg.config:3` hard-codes one user's `workDir`.
  - Pipeline-wide tuning (MINIMAP2_ALIGN memory doubling, RUN_SUMMARY_TABLES time) lives in
    the site profile.
  - The schema `$id`s are placeholder URLs.

### 3.7 Testing, tooling, repository hygiene

- **T1 · H · almost no automated tests.**
  - The only test, `tests/curate_models/test_curate_models.py`, passes when run directly,
    but it is a custom harness, one 263-line `main()`; pytest collects 0 tests.
  - No tests for 13 of 14 scripts (including `fomo_stats.py`, which owns the cross-script
    contracts) or for either core awk program.
  - No nf-test, no CI, no lint or format configuration, no pre-commit.

  Every "verified / byte-identical / fails loudly" claim in the docs was checked once, by hand.
- **T2 · H · the stub run is broken.** Every module has a stub, but `-stub-run -profile test`
  never completes:
  - without containers, version `eval`s fail (TD2);
  - with containers, the `CURATE_MODELS` stub runs `gzip`, which its image lacks;
  - the `GFF_BY_SOURCE` stub writes `allModels.gff3`, outside the `.projected.gff3` grammar,
    so the per-source path is never exercised and a misnamed file gets published.

  Unexercised code rots.
- **T3 · M · the test data cannot catch the documented production failure.** All committed
  test FASTAs are BGZF, while real assemblies are plain gzip. CLAUDE.md itself warns that a
  `faidx` shortcut "passes `-profile test` and fails on the cluster".
- **T4 · M · a vendored binary in git.** `bin/gff-feature-stats` is a 703 KB ELF, and each
  rebuild adds another blob. Its provenance (commit `67ffea7`; `-V` says 0.2.0 for two commits)
  lives only in prose. It forces the arm64 workaround and the unverified glibc floor.

  *Fix:* a tagged upstream release in a digest-pinned container.
- **T5 · L · hygiene.**
  - The tracked `.DS_Store`.
  - `.gitignore` entries for directories that do not exist.
  - A one-off dev script in `bin/` (staged into every task).
  - `plans/13_resource_rightsizing_audit.py`, a reusable tool, sits in `plans/`.
  - No CHANGELOG or CONTRIBUTING.
  - Large multi-plan commits with terse messages ("hide plans" ×2, "general results").

---

## 4. Principles

These rules address the root causes. Each comes with a mechanical check, so it does not
depend on discipline.

| # | principle | enforced by |
|:---|:---|:---|
| 1 | **Identity lives in meta, never in parsed strings.** Every (source, target, track) fact is a meta key set once. Filenames are *derived* from meta, never parsed back. A task that needs per-file identity reads a manifest. | grep-based CI check: no `.from_` slicing or `${meta.decoy ? …}` outside the helper |
| 2 | **One registry per vocabulary.** Tracks (`feature_type`, decoy, label, order, td2), MultiQC sections (suffix, id, order) and curation metric keys are each defined once (YAML/JSON) and read by Groovy and Python. | contract tests that round-trip every grammar and check every `sp` suffix has a producer |
| 3 | **`modules/nf-core` means unmodified and tool-managed.** Anything FOMO-specific is local and cut down to what is called. | `nf-core modules lint` in CI; `modules.json` present |
| 4 | **Invariants are tests, not prose.** A "must agree" or "fails loudly" statement either has a test or is deleted. | pytest fixtures per script, an nf-test stub matrix over (F, D, S, T), CI |
| 5 | **Comment policy: present tense, why, short.** At most ~3 lines unless it states an invariant. No run, job, clade or measurement anecdotes; no history; no file:line references; no `plans/NN` in user-facing text. History goes in commits, CHANGELOG and plans with a status. | `scripts/code_metrics.py` in CI with a budget (Nextflow layer ≤ 0.4, zero incident and history hits) |
| 6 | **One pattern per kind of thing.** One gating mechanism for optional features, one stats-producer shape (producer + adapter in one task), one batch helper, one scratch helper, three pinned images, one CLI convention. | review checklist in CONTRIBUTING |
| 7 | **Small, single-purpose commits.** One plan or change per commit; message = what + why. | PR template |

---

## 5. Roadmap

Each phase can be merged on its own, and each later phase relies on the safety net from
phase 0. Estimates are rough engineering days.

### Phase 0 — safety net (2–3 days)

**Do this before anything else.**

1. Convert `test_curate_models.py` to pytest (parametrized fixtures). Add fixture tests for
   `fomo_stats`, `gff_by_source`, `filter_gff_by_id`, `bam_to_gff.sh` and
   `filter_transcript.sh` on tiny inputs: multi-exon, trans-spliced, ties, `Parent=a,b`,
   empty input, duplicate IDs.
2. Fix the three stub defects (T2). Then add an **nf-test stub matrix** — lncRNA, +mRNA,
   +decoy, a pure target with and without gff3 — that asserts task counts against the
   formulas in CLAUDE.md and the published paths.
3. Add one plain-gzip FASTA and one uncompressed GFF3 to the test sheet (T3).
4. Add `pyproject.toml` (ruff lint and format; C901 ≤ 15), `shellcheck`, `nextflow lint`,
   `.editorconfig` and pre-commit. Add a GitHub Actions job running the above plus
   `-profile test,docker`.

### Phase 1 — delete and fix (1–2 days)

1. Everything in §2.
2. Fix the latent bugs: R1 (batch strict mode + retry codes), R2 (`curation_index`), R3
   (`bam_to_gff` source check), the species pattern (P1), and the two strict-parser errors
   (R4).
3. Rewrite user-facing text in present tense (schema help, MultiQC descriptions,
   `pipeline_info` naming — and either timestamp the traces or fix the comment).
4. One comment sweep applying principle 5: delete V1–V3, cap V4.

### Phase 2 — single sources of truth (3–4 days)

1. A track registry → `meta.track`, `meta.label`, `meta.order`, `meta.td2`, stamped in
   PREPROCESSING. Replace the 35 inline expressions and 5 ladders. Grow `fomo_stats.py` into
   the naming module; a Groovy helper covers what meta cannot.
2. Batch tasks write manifests; delete the three filename parsers (P1).
3. A MultiQC section registry plus its consistency test (P3).
4. Config: nf-core include order (P5), `'.*:NAME$'` selectors (P6), resources in one file with
   honest labels (C4), the site profile holding only site settings (R7).

### Phase 3 — consolidate modules (3–4 days)

1. Fork → local: `minimap2_align_splice`, `extract_transcripts` (C1). Reinstall the unchanged
   modules with `nf-core modules install` (C2).
2. GFF stats ×3 → one module; FILTER_ANNOTATION/FILTER_GFF_BY_ID → one; adapters run in the
   producing task (D4).
3. Shared `for_each_projected.sh` and scratch helper (D3, D4). Three pinned images (C5).
   Settle versions (C3).
4. Use meta joins and an id-keyed PREPROCESSING (P7). Drop the duplicate target stats (§2).

### Phase 4 — Python structure (3–5 days)

1. `fomo_gff.py` and `fomo_mqc.py` (D1, D2).
2. Split `curate_models.main` and `read_reference` (P8); dataclasses for Model, Gene and
   Transcript, plus a `ChainKey` NamedTuple.
3. `run_summary_tables`: cached `Inputs`, compute/emit per section, then data-JSON inputs
   instead of presentation TSVs (P4).
4. CLI and logging conventions (N4, R5); `relocate_loci` with a local RNG and incremental
   weights (R6).
5. `filter_transcript.sh` → Python on `fomo_gff` (P9).

### Phase 5 — documentation (1–2 days)

1. CLAUDE.md → ~150–200 lines: glossary, roles, present-tense invariants (each pointing to the
   code or test that enforces it), and three checklists — *add a track*, *add a stats section*,
   *add a param*.
2. README as the human entry point: requirements, quick start with an engine profile,
   samplesheet example, key params, outputs. Move usage, output and architecture detail into
   `docs/`.
3. `plans/`: an ADR header on each (Status: Accepted / Superseded by NN; Date; a 3-line
   decision) plus `plans/README.md` as an index. Add `CHANGELOG.md` with the "results not
   comparable before vX" flags. Add CONTRIBUTING with the principles in §4.

### What success looks like

**Adding a feature type** (e.g. `snoRNA`):

- **Today:** edits in ≥ 8 files across Groovy, config, bash and Python, and the MultiQC
  YAMLs. A missed copy drops rows silently.
- **After phases 0–2:** one registry entry, one fixture, and the stub-matrix test fails until
  every count and path is right.

**Adding a stats section:**

- **Today:** a producer, an adapter, a suffix, an `sp` pattern, an order key, a clean-ext and
  a `pick()` call, kept in sync by reading CLAUDE.md.
- **After:** a producer that writes data JSON, plus one registry line. The consistency test
  covers the rest.

---

## 6. Verification log

Re-verified by hand (on top of the reviewers' own checks):

| claim | check |
|:---|:---|
| track expression 35×, `meta.gtype` `_decoy` used once | `grep -rho` over `.nf`/`.config` |
| 3 ordering ladders in config, parsed back in Python | `modules.config:88,139,269`; `run_summary_tables.py:49,374,429,773` |
| `gunzip`, `gffcompare` unused | every `include` in workflows/subworkflows |
| strict-parser failure | `nextflow lint` (6 errors) and `NXF_SYNTAX_PARSER=v2 nextflow run … -preview` → "Script compilation failed" |
| batch masking, exit 123 | `bash`/`xargs` reproduction: `false; echo` → 0, `kill -9` → 123 |
| stub run broken (3 defects) | `-stub-run -profile test` with and without singularity (scratch dir) |
| `curation_index` precedence vs docstring | `run_summary_tables.py:214-222`, `curation_to_mqc.py:18-23` |
| byte-identical `attr()`, duplicated `TX_TYPES` | `diff`, `grep` |
| unused `is_aggregate`/`is_self`/`build_shape`, `--section-id` log-only | `grep` |
| stale `FILTER_TARGET`, `ANNOT_ORDER`, line refs, "schema rejects gff3", timestamps vs `overwrite` | the cited lines |
| re-keyed `id` unused | `reporting.nf:24-27`, `run_summary.nf:33-35` |
| include order | `nextflow.config:75-123` |
| no `topic('versions')` consumer, unused `big_mem`/`small` labels | `grep` |
| `bam_to_gff.sh` defaults source to `.` | `bam_to_gff.sh:77-80` |
| test FASTAs are BGZF; `samplesheet_full.csv` and 2 species unused; README test command | `file`, `grep` |
| custom test passes; pytest collects 0 | ran both |
| ruff complexity list | `ruff check --select C901,PLR0912,PLR0915` |

Reported by a reviewer and plausible but not re-run by me:

- `filter_transcript.sh` ignores exons written as `; Parent=` (with a space).
- The mawk `--version` output.
- The upstream nf-core diffs.
- `nextflow config` confirming the include-order fix.

Tools: ruff 0.11.5, pyflakes 3.4, Nextflow 25.10.4 (`lint`, v2 parser); metrics from
`scripts/code_metrics.py` (`python3 ablation/scripts/code_metrics.py` from the repo root).
