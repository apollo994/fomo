# 18 — Per-target batching: one process for all sources

## Status: implemented (2026-09-08), verified against baseline on `-profile test`

Implemented as `modules/local/gff_stats_projected_batch.nf` (`GFF_STATS_PROJECTED_BATCH`)
and `modules/local/gffcompare_batch.nf` (`GFFCOMPARE_BATCH`), wired into
`subworkflows/local/projection.nf` and `subworkflows/local/benchmarking.nf`. Two things in
the section below turned out to be **wrong as originally sketched** and were fixed during
implementation — both caught by actually diffing a real run against baseline, not by
inspection:

1. **"The naming grammar already round-trips" was only true for the S per-source split
   files, not the 2 aggregates.** `allModels_raw`/`allModels_collapsed`'s physical
   filenames (`<target>.allModels.<gtype>.raw|collapsed.gff3`, from
   `FILTER_ALLMODELS`/`COMBINED_GTF_TO_GFF`) never followed the `.from_<id>.` convention —
   the OLD per-source `ext.prefix` computed their output name purely from `meta.id`,
   independent of that input filename, so `${gff%.gff3}` doesn't recover it. Fix: both new
   modules symlink the 2 aggregate files to the same `.from_allModels_raw|collapsed.<gtype>.
   projected.gff3` shape as the per-source files (using the two known pseudo-source ids
   directly — `meta` is available in the script), and glob only `*.projected.gff3` for the
   main loop.
2. **`GFF_STATS_PROJECTED`'s old prefix included `.projected`; `GFFCOMPARE`'s did not.**
   So `GFF_STATS_PROJECTED_BATCH` strips `.gff3` (prefix keeps `...projected`), while
   `GFFCOMPARE_BATCH` strips the longer `.projected.gff3` (prefix drops it) — verified by
   diffing `.gffcompare.stats` filenames against baseline; getting this wrong broke
   `fomo_stats.py`'s `FNAME_RE` silently (files parsed as unrecognised, dropped from
   ranking) rather than erroring.
3. **The old `GFF_STATS_PROJECTED_TO_MQC`'s `ext.sample_name` ordering token
   (`1_lncRNA`/`2_decoy_lncRNA`/…) was dropped in the first pass**, since the fused
   module derived `--name` from the file prefix directly. `bin/run_summary_tables.py`
   parses that exact token back out (`ANNOT_ORDER` map) to key its per-target/per-class
   aggregates — losing it didn't error, it silently nulled out most of
   `run_summary.json` and left `Lycaena_hippothoe_580924: {}` empty. Fix: `order` is
   computed once per task (constant — depends only on `feature_type`/`decoy`, both
   task-level), the per-file `id` is recovered from the `.from_<id>.` filename segment
   inside the loop (`${gff#*.from_}` / `${id%%.*}`), and `--name` is rebuilt as
   `${target_id}.from_${id}.${order}` — byte-identical to the old sample names.

**Verification performed:** `nextflow lint .` clean on every touched file; a real
`-profile test,singularity` run (not `-stub-run`) completed end-to-end (76 tasks vs. 93
baseline); diffed in full against a baseline run from `main` (via `git stash`). After fix
#3, `run_summary.json`, every `*_mqc.tsv`/`.json` table, `select_top_sources.py`'s ranking
output, and every `.gffcompare.stats`/`.stats.json` numeric payload are **byte-identical**
to baseline. The only remaining, accepted differences: the 2 aggregate items' `.tmap`/
`.refmap` filenames (which embed the query file's own name — now the symlink name, not
the original) and the corresponding `.gffcompare.stats` file's self-documenting `#gffcompare
...`/`#= Summary for dataset:` header lines (which echo that same filename) — cosmetic,
verified to carry no numeric difference, and invisible to every downstream parser (which
reads the *filename*, or specific data lines, never that header).

## Problem

The Lepidoptera all-vs-all run (S=531, T=1260, T_g=531, `-resume 1f752880-…`) has OOM-killed
the Nextflow **head** JVM four times in a row on the interactive `srun` shell it runs in,
despite requested memory escalating 2 GB → 2 GB → 8 GB → 8 GB → **16 GB** (`sacct`, all four
non-trivial attempts show `State=OUT_OF_MEMORY`). The most recent attempt (job `28291750`,
16 GB) reached `MaxRSS = 16.38 GB` — effectively the cgroup ceiling — in under 6 hours,
against 16d 20h for the previous 8 GB attempt: the trend is not "give it more RAM," it is
that the graph itself is too large for one JVM to track through a `-resume`.

Per `CLAUDE.md`'s own cardinality accounting, two stages account for 952k of the run's ~968k
tasks:

| Stage | Formula | At S=531, T=1260, T_g=531 |
|---|---|---|
| `PROJECTION:GFF_STATS_PROJECTED` (+ `_TO_MQC` adapter) | `2·F·D·S·T` | **1,338,120** (measured partial: 41,376 + 19,585 submitted before the last kill) |
| `BENCHMARKING:GFFCOMPARE` | `(F·S·D + 2·F·D)·T_g` | **283,023** |

`plans/17_scheduling_efficiency.md` explicitly named this "the quadratic task counts that
make a 1260-species run infeasible regardless of latency" and scoped it out. This plan is
that follow-up.

Per-task cost is trivial — one `gff-feature-stats` call is ~1 s at ~20 MB RSS
(`modules/local/gff_stats.nf:5-8`), one `GFFCOMPARE` is a few seconds at <1 GB
(`conf/base.config:25`). The problem is **granularity**: ~1.6M SLURM submissions and,
worse, ~1.6M `TaskHandler` objects the head JVM must carry through a `-resume`, plus their
proportional contribution to the `mqc_files`/`stats` unions that `RUN_SUMMARY` and
`CONSENSUS_TOP` `.collect()` in full before their own single downstream task — very likely
why the 16 GB attempt spiked so fast if this resume's DAG frontier had already reached those
collectors.

## What changes

**One task per target does every source**, instead of one task per (source, target) pair.
Both `SPLIT_GFF_BY_SOURCE` (already: one task emits all S files for a target — see
`subworkflows/local/projection.nf:170-172`) and the redesigned stats/gffcompare stages loop
internally over every source's file with `xargs -P ${task.cpus}`, so the same wall-clock work
that used to be S separate SLURM jobs becomes one job using S/`task.cpus` sequential rounds.
This is embarrassingly parallel — no source's `gff-feature-stats`/`gffcompare` call touches
another's — so cores map directly to wall-clock, not just to throughput.

### Why not gffcompare's native multi-query mode

`gffcompare -r ref.gff3 q1.gff3 q2.gff3 …` is tempting — the nf-core module's `path(gtfs)`
input already accepts a list (`modules/nf-core/gffcompare/main.nf:11,39`) — but it **pools**
every query's transcripts into one shared comparison, reporting one aggregate Sn/Pr across all
of them. `CONSENSUS_TOP:SELECT_TOP_SOURCES` needs a **per-source** F1 to rank donors
(`subworkflows/local/consensus_top.nf:22-49`); a pooled comparison would destroy exactly the
signal the ranking depends on. So batching here means **looping S+2 separate single-query
invocations inside one task**, not merging them into one gffcompare call. Each loop iteration
still writes its own complete `.stats`/`.tracking`/`.loci`/`.tmap`/`.refmap` under its own
`-o` prefix — nothing about gffcompare's per-file output changes, only how many times it's
invoked per SLURM job.

### The naming grammar already round-trips — no manifest needed

Every consumer of a batch item needs to know *which* output name to write. It already does,
for free: `ext.prefix` for both processes today is built as
`"${meta.target_id}.from_${meta.id}.${feature_type}${decoy}"` — and that string is *exactly*
the basename `SPLIT_GFF_BY_SOURCE`'s `--name-template` already gave the input file
(`<target>.from_{source}.<gtype>.projected.gff3`, `projection.nf:170-172`, filename grammar
documented in `bin/gff_by_source.py:17`). The two aggregate models
(`allModels_raw`/`allModels_collapsed`) go through the identical `.from_<id>.` construction
because `meta.id` is just set to those literal strings before this stage
(`projection.nf:222-223`). So each loop iteration can derive its own output prefix as
`${gff%.gff3}` — **no per-item manifest, no Groovy meta needed inside the script.**

### Glob trap to avoid

The batch's work dir also holds the target's one reference GFF3 (for `GFFCOMPARE_BATCH`),
whose `ext.prefix` is `${meta.id}.${meta.feature_type}` (`conf/modules.config:300`) — it does
**not** end in `.projected.gff3` / `.raw.gff3` / `.collapsed.gff3`, but a bare `*.gff3` glob
would still catch it if anyone loosens the pattern later. Loop over the three known query
suffixes explicitly, never a bare `*.gff3`:
```bash
printf '%s\0' *.projected.gff3 *.raw.gff3 *.collapsed.gff3 | xargs -0 -P "${task.cpus}" ...
```

### Cardinality

| Stage | Before | After |
|---|---|---|
| `GFF_STATS_PROJECTED` + `_TO_MQC` (fused, see below) | `2·F·D·S·T` ≈ 1,338,120 | **`F·D·T`** ≈ 1,260 |
| `BENCHMARKING:GFFCOMPARE` | `(F·S·D+2·F·D)·T_g` ≈ 283,023 | **`F·D·T_g`** ≈ 531 |

Both collapse by a factor of ~S (≈531×). Total task count for the run falls from ~968k to
roughly **~15k** (dominated now by `PREPROCESSING` and the alignment chain, which plan 15
already made linear in T).

### Decisions

| Question | Decision | Rationale |
|---|---|---|
| Batch scope | **All S sources + the 2 aggregates, per (target, feature_type, decoy)** — not a tunable batch size | Matches the ask directly; `SPLIT_GFF_BY_SOURCE` already materialises the full per-target file set in one task, so there is no natural sub-batch boundary to pick, and picking one just adds a knob with no benefit here |
| Intra-task parallelism | `xargs -P ${task.cpus}`, not GNU `parallel` | Always present (coreutils); avoids adding a package to two containers that don't ship it today |
| CPU/mem for the two new processes | Reuse the existing, currently-unused `process_high` label (`conf/base.config:42-47`: 12 cpu / 72 GB / 2 h) as the starting point, right-size later | No new config block; matches this repo's own "trim later, from real numbers" convention (`plans/13_resource_rightsizing.md`) |
| Failure handling | Fail-fast — `xargs` (no `-k`) aborts the batch task on the first non-zero exit, Nextflow retries/fails the whole target's stage | Matches `bin/gff_by_source.py`'s own philosophy ("fails loudly rather than corrupting"); per-source failures should be structurally impossible here since `gff_by_source.py` already validated the input upstream |
| Fuse `GFF_STATS_PROJECTED` + `_TO_MQC`? | **Yes** — one task calls `gff-feature-stats` then `gff_stats_to_mqc.py` per source, JSON never leaves the task | Removes the adapter's task count entirely rather than merely batching it; the JSON was never published (not in `CLAUDE.md`'s "Result publishing" list), so nothing external changes. `GFF_STATS_TARGET`/`GFF_STATS_TOP` (small counts, `F·T_g` / `2·F·T_g`) are **not** touched — the two-step convention stays for everything that isn't S·T-scaled |
| Fuse `GFFCOMPARE`'s outputs? | No — keep the existing 5 output globs (`.stats`/`.tracking`/`.loci`/`.tmap`/`.refmap`), all of which are published (`conf/modules.config` "Result publishing" section) | Nothing to fuse; each loop iteration already writes its own complete set under its own `-o` prefix, unchanged from today |

## Downstream contract: what has to keep working

Traced every consumer of the two changed emits before touching anything:

- `PROJECTION.out.gff3` is consumed **only** by `BENCHMARKING` (`workflows/fomo.nf:77`) — safe
  to reshape from "S+2 individual items" to "one per-target list", since I'm rewriting both
  ends together.
- `BENCHMARKING.out.stats` is consumed by `CONSENSUS_TOP` (`ch_all_stats`,
  `consensus_top.nf:17`) and by `REPORTING`/`RUN_SUMMARY` via the `mqc_files` union.
  `CONSENSUS_TOP:SELECT_TOP_SOURCES`'s ranking filter
  (`consensus_top.nf:40-45`) needs **per-source** `meta.id`/`meta.decoy`/`meta.feature_type`
  on each stats file — to exclude the two aggregate-id families and self
  (`meta.id != meta.target_id`) **before** grouping per target. So `GFFCOMPARE_BATCH`'s
  output (one shared target-level meta, many files) must be exploded back into individual
  `(meta, statsFile)` pairs with the *correct* re-derived `id`, using the exact same
  string-slice re-key already proven for the GFF3 split
  (`projection.nf:178-187`):
  ```groovy
  ch_stats_persource = GFFCOMPARE_BATCH.out.stats
      .transpose()
      .map { meta, stats ->
          def pre = "${meta.target_id}.from_"
          def suf = ".${meta.feature_type}${meta.decoy ? '.decoy' : ''}.gffcompare.stats"
          tuple([id:           stats.name[pre.size()..-(suf.size() + 1)],
                 feature_type: meta.feature_type,
                 decoy:        meta.decoy,
                 target_id:    meta.target_id], stats)
      }
  ```
  This is the one place the batching is **not** transparent downstream — everything else
  (MultiQC `sp` patterns, `bin/fomo_stats.py`'s filename parser) already routes by filename,
  not by Nextflow task/channel-item identity, so it doesn't care how many files came out of
  one task.
- `GFF_STATS_PROJECTED_BATCH`'s `mqc_files` do **not** need this re-key: `REPORTING` routes
  solely on `meta.target_id` presence (`CLAUDE.md`, "One report per target"), which is
  identical for every file in a batch, so a single shared meta is sufficient.
- `bin/fomo_stats.py`'s `AGGREGATE_IDS`/`FNAME_RE` (already the ground truth for filename
  parsing) do not change — the physical filenames are unchanged; only which task produced
  them changes. **This plan adds a third place that must agree with those two** (the two
  already named in `CLAUDE.md`: the channel filter in `consensus_top.nf`, and
  `AGGREGATE_IDS` in `fomo_stats.py`) — the re-key closure above. Comment it as such at the
  call site.

## Implementation phases

### Phase A — `GFF_STATS_PROJECTED` batching + fusion (lower risk: no re-key needed)

1. **New module** `modules/local/gff_stats_projected_batch.nf`, process
   `GFF_STATS_PROJECTED_BATCH`, `label 'process_high'`. Container: switch from
   `ubuntu:22.04` to `python:3.11` (Debian bookworm base) so both the vendored
   `gff-feature-stats` binary (needs glibc ≥ 2.34) and `gff_stats_to_mqc.py` (needs
   python3) are satisfied by one image — **verify** `python:3.11`'s glibc with the same
   smoke test `CLAUDE.md` already prescribes for the Rust binary
   (`docker run python:3.11 ldd --version`, or the `singularity exec ... -V` pattern).
   Input `tuple(val(meta), path(gffs))` where `gffs` is the whole per-target list (every
   source + the 2 aggregates for one `(target_id, feature_type, decoy)`). Output
   `tuple(val(meta), path("*_gffstats_projection_mqc.tsv"), emit: tsv)` plus the usual
   `versions` topic (unchanged content — version doesn't vary per source, so one emit per
   task is *more* correct than before, not less). Script sketch:
   ```bash
   run_one() {
       gff="$1"; prefix="${gff%.gff3}"
       gff-feature-stats "$gff" "${prefix}.stats.json"
       gff_stats_to_mqc.py --input "${prefix}.stats.json" --name "$prefix" \
           --role projection --kind projected --section-id gffstats_projection \
           --output "${prefix}_gffstats_projection_mqc.tsv"
   }
   export -f run_one
   printf '%s\0' *.projected.gff3 *.raw.gff3 *.collapsed.gff3 \
       | xargs -0 -P "${task.cpus}" -I{} bash -c 'run_one "$@"' _ {}
   ```
2. **`subworkflows/local/projection.nf`**: replace the `GFF_STATS_PROJECTED` +
   `GFF_STATS_PROJECTED_TO_MQC` includes/calls (`:8-9`, `:226-229`) with one
   `GFF_STATS_PROJECTED_BATCH` call. Feed it the **per-target list**, built directly from
   `SPLIT_GFF_BY_SOURCE.out.gff3` (already a List per `(target, feature_type, decoy)` before
   `.transpose()`) plus the two aggregate GFF3s folded in as extra list members, keyed on
   `[target_id, feature_type, decoy]`. Do **not** call `.transpose()` for this consumer at
   all — that's the entire point.
3. `PROJECTION.emit.gff3` — **decide alongside Phase B** whether to keep emitting the
   per-source `.transpose()`d form (still cheap here, ~S·T channel items but zero extra
   *tasks*) purely for BENCHMARKING to re-batch itself, or emit the per-target list directly.
   Emitting the list directly is the cleaner interface and is what Phase B assumes below —
   do Phase A and B together if convenient, or land Phase A first with the emit unchanged
   (keep `.transpose()` just for the emit, not for `GFF_STATS_PROJECTED_BATCH`) and do the
   emit-shape change in Phase B. Record which was chosen.
4. **`nextflow.config`** (docker profile arm64 pin): the selector
   `'.*:GFF_STATS(_TARGET|_PROJECTED|_TOP)?$'` is `$`-anchored and will **stop matching**
   `GFF_STATS_PROJECTED_BATCH` — update it to
   `'.*:GFF_STATS(_TARGET|_PROJECTED(_BATCH)?|_TOP)?$'`. This is exactly the kind of
   silent-drift trap `CLAUDE.md` calls out for this selector; do not skip it.
5. **`conf/modules.config`**: delete the now-dead `FOMO:PROJECTION:GFF_STATS_PROJECTED`/
   `_TO_MQC` selectors (`:273-281`), add one for `GFF_STATS_PROJECTED_BATCH` (no `ext.prefix`
   needed — the script derives names itself; keep `ext.role`/`ext.section_id` equivalents
   as hardcoded script constants or pass through `ext.args` if a future track needs to vary
   them).
6. **`CLAUDE.md`**: update the cardinality line for `GFF_STATS_PROJECTED` and add a short
   note next to the existing `RUN_SUMMARY`/`mqc_files` prose about the batching.

**Checkpoint A**: `-profile test,docker` run completes; `*_gffstats_projection_mqc.tsv`
content byte-identical to the pre-change run (diff `results/` against
`results_baseline/`); `grep -c GFF_STATS_PROJECTED_BATCH execution_trace.txt` equals
`F·D·T` for the test samplesheet, not `S·T`.

### Phase B — `BENCHMARKING:GFFCOMPARE` batching + re-key

1. **New module** `modules/local/gffcompare_batch.nf`, process `GFFCOMPARE_BATCH`,
   `label 'process_high'`, container reused verbatim from the nf-core module
   (`quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0` /
   `depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9f5acd7_0` — copy the
   `containerEngine`-conditional block). Inputs: `tuple(val(meta), path(gffs))` (the
   per-target query list, same shape Phase A already built) and
   `tuple(val(ref_meta), path(reference))` (the target's one filtered reference — a plain
   1:1 value now, not a `combine(by:...)` fan-out). Outputs: the same five globs the
   nf-core module already declares, unchanged (`*.stats`, `*.tracking`, `*.loci`, `*.tmap`
   optional, `*.refmap` optional, `*.annotated.gtf` optional). Script sketch:
   ```bash
   run_one() {
       gff="$1"; prefix="${gff%.gff3}"
       gffcompare -M --no-merge -r "${reference}" -o "${prefix}.gffcompare" "$gff"
       mv "${prefix}.gffcompare" "${prefix}.gffcompare.stats"
   }
   export -f run_one
   printf '%s\0' *.projected.gff3 *.raw.gff3 *.collapsed.gff3 \
       | xargs -0 -P "${task.cpus}" -I{} bash -c 'run_one "$@"' _ {}
   ```
2. **`subworkflows/local/benchmarking.nf`**: replace the per-source `combine(by:[0,1])` +
   `multiMap` + single-query `GFFCOMPARE` (`:55-72`) with a 1:1 join on
   `[target_id, feature_type]` between the per-target query list (from `PROJECTION.out.gff3`,
   now list-shaped per Phase A step 3) and `FILTER_TARGET.out.gff3` — no fan-out needed,
   since there is exactly one reference and one query-list per target.
3. Add the **re-key transpose** described above (`ch_stats_persource`) right after
   `GFFCOMPARE_BATCH`, and emit *that* as `BENCHMARKING.out.stats` — this is the piece
   `CONSENSUS_TOP` depends on; do not emit the raw batched output directly.
4. **`conf/modules.config`**: delete the dead `FOMO:BENCHMARKING:GFFCOMPARE` selector
   (`:303-306`), add one for `GFFCOMPARE_BATCH` carrying the same `ext.args = '-M --no-merge'`
   (or hardcode it in the script — either is fine since it no longer varies per item).
5. **`CLAUDE.md`**: update the `GFFCOMPARE` cardinality line and the "Multi-target
   wiring" note about join keys (still keyed on `[target_id, feature_type]`, just 1:1 now
   instead of fanned out).

**Checkpoint B**: same `-profile test,docker` diff as Checkpoint A, plus explicitly diff
`select_top_sources` output (`<target>.top_sources.csv`) byte-for-byte against baseline —
this is the one output whose *correctness*, not just task count, depends on the re-key step
picking up the right `meta.id` per file. A wrong re-key would silently rank sources under
the wrong name or leak an aggregate into the pool; the CSV diff is the only place that would
show up.

## Files touched

| File | Change |
|---|---|
| `modules/local/gff_stats_projected_batch.nf` | new — fused, batched, `xargs -P` |
| `modules/local/gffcompare_batch.nf` | new — batched, `xargs -P`, same 5 output globs as the nf-core module |
| `modules/local/gff_stats.nf`, `modules/local/gff_stats_to_mqc.nf` | unchanged — still used by `GFF_STATS_TARGET`/`GFF_STATS_TOP` (small, not S·T-scaled) |
| `modules/nf-core/gffcompare/main.nf` | unchanged — still used by `GFFCOMPARE_COMBINE`, `GFFCOMPARE_COMBINE_TOP`, `GFFCOMPARE_TOP` (single-input or already-small) |
| `subworkflows/local/projection.nf` | drop `GFF_STATS_PROJECTED`/`_TO_MQC` includes+calls; add `GFF_STATS_PROJECTED_BATCH`; stop transposing before stats; reshape `gff3` emit to per-target lists |
| `subworkflows/local/benchmarking.nf` | drop `GFFCOMPARE` include+call; add `GFFCOMPARE_BATCH`; 1:1 join instead of fan-out `combine`; add the re-key transpose |
| `conf/modules.config` | delete 2 dead selectors, add 2 new ones (no `ext.prefix` needed on either) |
| `nextflow.config` | widen the arm64 `GFF_STATS*` docker-platform-pin selector to match `_BATCH` |
| `CLAUDE.md` | cardinality table, subworkflow map emits, "third place to keep in sync" note next to the existing `AGGREGATE_IDS` callout |

## Verification

No nf-test harness exists in this repo yet (per `plans/17`'s own verification section); the
profile smoke test is the regression net.

1. **Baseline first.** `nextflow run . -profile test,docker` on current `main`, copy
   `results/` → `results_baseline/`.
2. **After each phase**, re-run the same command and diff. Every published filename and file
   *content* must be byte-identical to baseline — the naming grammar is untouched by
   construction (§ "already round-trips" above), so any diff means a real bug, not a
   deliberate rename.
3. **Task-count assertion** on `results/pipeline_info/execution_trace.txt`:
   `grep -c GFF_STATS_PROJECTED_BATCH` and `grep -c GFFCOMPARE_BATCH` must equal `F·D·T` and
   `F·D·T_g` respectively for the test samplesheet (5 species) — i.e. exactly one task per
   target per (feature_type, decoy) combination, not one per species.
4. **`select_top_sources` CSV diff** (Checkpoint B) — the correctness-sensitive check for the
   re-key step.
5. **Real-cluster validation before touching the live 1260-species run**: 3–5 rows of the
   real Lepidoptera samplesheet (at least one `both`, one `target`) under `-profile crg`,
   watching `sacct` for the new processes' actual `MaxRSS`/`Elapsed` to right-size
   `process_high`'s 12 cpu / 72 GB down from the placeholder (per `plans/13`'s
   audit-then-trim convention) — 12 cpus × ~20 MB RSS for `gff-feature-stats` or a few
   hundred MB for `gffcompare` should need nowhere near 72 GB; measure, don't guess.
6. **Head-memory win, the actual point of this plan**: kick off a `-resume` of the real
   1f752880 run afterward with a modest, fixed head allocation (e.g. 8 GB, `NXF_OPTS` with
   `-Xmx6g -XX:+HeapDumpOnOutOfMemoryError` per the earlier diagnostic recommendation) and
   confirm it no longer climbs toward the ceiling within hours. Everything upstream of the
   two changed stages stays cached under `-resume` — expect only `GFF_STATS_PROJECTED*` and
   `GFFCOMPARE`/`BENCHMARKING` (and everything transitively downstream: `CONSENSUS_TOP`,
   `REPORTING`, `RUN_SUMMARY`) to re-execute, since the task signature changes.

## Baseline

`results_baseline/` from `nextflow run . -profile test,docker` on current `main`, taken
before Phase A starts — both phases diff against it per the Verification section.
