# Plan: Resource right-sizing from execution-trace evidence

## Context

FOMO's resource tiers in `conf/base.config` were set before any profiling existed. Two
completed Lycaenidae runs were audited from their Nextflow traces:

- `…/fomo_runs/lycaenidae/targets/Polyommatus_icarus_265386/pipeline_info/execution_trace.txt` (1424 tasks, all `COMPLETED`, first attempt)
- `…/fomo_runs/lycaenidae/targets/Aricia_agestis_91739/pipeline_info/execution_trace.txt` (1424 tasks, 684 `COMPLETED` + 740 `CACHED`, first attempt)

Both agree closely, and no task in either run failed or retried. Measured against what was
requested:

| Resource | Reserved | Actually used | Overestimate |
|---|---|---|---|
| Memory | 21.5 / 24.0 GB·h | 3.07 / 3.36 GB·h | **7×** |
| CPU | 3.60 / 4.01 cpu·h | 1.46 / 1.31 cpu·h | **2.5–3.1×** |
| Walltime | 10 m / 1 h limits | median task **0.5 s**, longest **89 s** | 40–1000× |

Per-process worst offenders (max observed peak RSS vs. request, both runs agree within ~5%):

| Process | Label | n | Requested | Max peak RSS | Max cores | Max runtime |
|---|---|---|---|---|---|---|
| `PROJECTION:MINIMAP2_ALIGN` | `process_medium` | 80 | 24 GB / 4 cpu / 1 h | 6.8 GB | 3.66 | 89 s |
| `*:TD2_PREDICT` (both instances) | `process_medium` | 40 | 24 GB / 4 cpu / 1 h | 0.60 GB | 1.24 | 15 s |
| `PROJECTION:MINIMAP2_INDEX` | (`withName` in crg.config) | 1 | 32 GB / 8 cpu / 1 h | 6.7 GB | 1.54 | 22 s |
| `EXTRACT_SEQUENCES` / `EXTRACT_DECOY_SEQUENCES` / `EXTRACT_PROJECTED_LNC` / `COMBINED_GTF_TO_GFF` (all `GFFREAD`) | `process_low` | ~105 | 12 GB / 2 cpu / 10 m | 0.23 GB | 1.00 | 33 s |
| `PREPROCESSING:SEQKIT_STATS` | `process_low` | 80 | 12 GB / 2 cpu / 10 m | 0.03 GB | 0.40 | 2 s |
| everything else (~900 tasks: `GFF_STATS*`, `*_TO_MQC`, `BAM_TO_GFF`, `SAMTOOLS_*`, `RELOCATE_LOCI`, `GFFCOMPARE*`, `FILTER_*`, `GUNZIP*`, `RENAME_FASTA_HEADERS`, `BEDTOOLS_COMPLEMENT`, `MULTIQC`, `SELECT_TOP_SOURCES`) | `process_single` | ~900 | 6 GB / 1 cpu / 10 m | **0.78 GB** (worst: `GFFCOMPARE_COMBINE`) | 0.93 | 45 s |

Full audit methodology and a reusable script are in `plans/13_resource_rightsizing_audit.py`.

## Constraints the implementer MUST respect

Read this section before editing anything. Three of these will bite you.

1. **Retry escalation is multiplicative and already works.** `base.config` uses
   `memory = { N.GB * task.attempt }` with `maxRetries = 3` and
   `errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }`.
   SLURM OOM-kill (137) and time-limit kill (SIGTERM 143 / SIGKILL 137) both fall in that
   range, so an under-provisioned task retries with 2×/3×/4× the resources rather than
   failing the run. **Preserve the `* task.attempt` form on every value you change.** This
   is what makes tightening these tiers a recoverable change instead of a risky one.

2. **Changing `cpus` invalidates the Nextflow cache; changing `memory`/`time` does not.**
   Four modules interpolate `task.cpus` into the command line, so altering their CPU count
   changes the script text and therefore the task hash:
   - `modules/nf-core/minimap2/align/main.nf`
   - `modules/nf-core/minimap2/index/main.nf`
   - `modules/nf-core/samtools/stats/main.nf`
   - `modules/nf-core/seqkit/stats/main.nf`

   This matters more than usual here: `conf/crg.config` sets a **shared**
   `workDir = '/nfs/scratch01/rg/fzanarello/work_fomo'` with `cache = 'lenient'`, used by
   all ~20 Lycaenidae target runs. Invalidating `MINIMAP2_INDEX` in particular cascades to
   all 80 `MINIMAP2_ALIGN` tasks per target and everything downstream of them. The steps
   below are ordered so that steps 1–3 cause **zero** cache churn.

3. **Do NOT globally shrink minimap2 memory.** `conf/crg.config` documents, in a comment
   on the `MINIMAP2_INDEX` override, that indexing a ~3 GB genome with `-x splice` needs
   ~25–30 GB and is OOM-killed at 12 GB. Lycaenidae assemblies are ~0.4–0.6 Gb, which is
   why peak RSS here is only 5.9–6.8 GB. **Minimap2 memory scales with target genome size,
   so the 24/32 GB requests are not wrong in general — they are wrong for this clade.**
   The fix must be clade-scoped (step 4) or size-aware (step 6), never a global constant.

4. `MINIMAP2_INDEX` already passes `-t $task.cpus` (`index/main.nf:24`). The 1.15–1.54
   cores observed is minimap2 index construction being largely serial. There is no missing
   flag to fix — just an over-request to trim.

5. `process_high` (12 cpu / 72 GB / 2 h) is **unused** — no process in either trace
   requested it. Leave it alone; do not "clean it up".

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| `process_single` memory | 6 GB → **2 GB** | worst observed 0.78 GB → 2.6× headroom, 8 GB after retries |
| `process_low` memory | 12 GB → **4 GB** | worst observed 0.23 GB |
| `process_low` cpus | 2 → **1** | worst observed 1.00 core; `SEQKIT_STATS` pinned back to 2 to avoid cache churn |
| `process_medium` cpus | **unchanged (4)** | `MINIMAP2_ALIGN` genuinely uses 3.4–3.66 cores; changing it invalidates cache for no gain |
| `process_medium` memory | **unchanged (24 GB)** | genome-size dependent — handled clade-locally in step 4, generally in step 6 |
| `process_medium` time | 1 h → **20 min** | longest observed 89 s; improves SLURM backfill priority |
| `TD2_PREDICT` label | `process_medium` → **`process_low`** | 0.60 GB / 1.24 cores / 15 s — nothing medium about it |
| `MINIMAP2_INDEX` cpus | 8 → 2, **deferred to step 5** | correct but cache-invalidating; bundle with a planned full re-run |
| `process_single` time | **unchanged (10 min)** | already the sane floor for a SLURM job; ratio looks bad, absolute cost is nil |
| `process_high` | **unchanged** | unused by this pipeline |

Expected outcome per target: reserved memory **~22 GB·h → ~9 GB·h (−58%)**, reserved CPU
**~3.6 → ~2.6 cpu·h**, with steps 1–4 applied.

---

## 1. Retune `process_single` and `process_low` in `conf/base.config`

Edit the two `withLabel` blocks. Keep the `* task.attempt` multipliers exactly as shown.

```groovy
    withLabel: process_single {
        cpus   = { 1                       }
        memory = { 2.GB       * task.attempt }   // was 6.GB  — worst observed peak 0.78 GB (GFFCOMPARE_COMBINE)
        time   = { 10.minutes * task.attempt }   // unchanged
    }
    withLabel: process_low {
        cpus   = { 1                       }     // was 2 * task.attempt — worst observed 1.00 core
        memory = { 4.GB       * task.attempt }   // was 12.GB — worst observed peak 0.23 GB (GFFREAD)
        time   = { 15.minutes * task.attempt }   // was 10.minutes — GFFREAD hit 33 s, headroom for larger clades
    }
```

Note the `cpus` closures deliberately drop `* task.attempt`: escalating CPUs on retry does
nothing for an OOM or a timeout on a single-threaded tool, and for the four `task.cpus`
modules it would produce a different command on every attempt.

Leave the top-level default block (`cpus = 1`, `memory = 6.GB`, `time = 10.minutes`)
untouched as a fallback for any future unlabelled process.

## 2. Pin `SEQKIT_STATS` back to 2 CPUs to avoid cache invalidation

`SEQKIT_STATS` is `process_low` and interpolates `task.cpus` into its command, so step 1's
`cpus 2 → 1` would re-run all 80 instances per target plus `SEQKIT_TO_MQC` and `MULTIQC`
downstream. It only uses 0.40 cores, but the saving (0.03 cpu·h) is not worth the churn.
Add to the `process` block in `conf/base.config`, after the `withLabel` blocks:

```groovy
    // SEQKIT_STATS interpolates task.cpus into its command line; holding cpus at 2 keeps
    // the task hash stable against the shared workDir. It only uses ~0.4 cores — revisit
    // (together with MINIMAP2_INDEX, see plans/13) at the next intentional full re-run.
    withName: '.*:SEQKIT_STATS' {
        cpus = 2
    }
```

Verify after implementing that this selector actually matches — the trace shows the fully
qualified name `FOMO:PREPROCESSING:SEQKIT_STATS`. Note `conf/modules.config` already uses
the unanchored form `withName: 'FOMO:PREPROCESSING:SEQKIT_STATS'`; either form is fine, but
be aware that file's comment at line ~54 documents that `withName` matching needs `$`
anchors when one process name is a prefix of another (`GFF_STATS` vs `GFF_STATS_TO_MQC`).
`SEQKIT_STATS` has no such collision.

## 3. Move `TD2_PREDICT` off `process_medium`

In `modules/local/td2_predict.nf`, line 3:

```groovy
-    label 'process_medium'
+    label 'process_low'
```

Both instances (`PREPROCESSING:TD2_NONCODING:TD2_PREDICT` and
`PROJECTION:TD2_NONCODING_PROJ:TD2_PREDICT`) peaked at 0.60 GB / 1.24 cores / 15 s against
a 24 GB / 4 cpu / 1 h request. A label change alters no command text, so there is no cache
impact.

Caveat to record in the commit message: PSAURON/TD2 memory scales with transcript count,
and Lycaenidae lncRNA sets are small. With `process_low` at 4 GB and multiplicative retry,
a larger clade gets 4 → 8 → 12 → 16 GB before failing. If a future clade retries here
routinely, give `TD2_PREDICT` its own `withName` block rather than reverting the label.

Also leave the commented-out GPU `withName: '.*:TD2_PREDICT'` block in `conf/crg.config`
exactly as-is — it is a documented future path, and a GPU build would change the CPU/memory
picture entirely.

## 4. Cap minimap2 memory for this clade only (run-level config, no repo change)

This is where 63% of the wasted memory·hours sit, and it is the one number that must not be
hard-coded globally (see Constraint 3). Create a new file in the run directory:

`/users/project/phyloENCODE_009142_no_backup/general_analyses/crossmapping/fomo_runs/lycaenidae/resources_lycaenidae.config`

```groovy
// Clade-scoped minimap2 memory for Lycaenidae (~0.4-0.6 Gb assemblies).
// Measured peak RSS: MINIMAP2_ALIGN 6.8 GB, MINIMAP2_INDEX 6.7 GB across 2 audited targets
// (see fomo/plans/13_resource_rightsizing.md). 12 GB leaves ~1.8x headroom; retry
// escalation reaches 48 GB if a larger assembly sneaks into the samplesheet.
// DO NOT promote these values into conf/base.config or conf/crg.config — minimap2 memory
// scales with target genome size and a ~3 Gb genome needs 25-30 GB.
process {
    withName: 'FOMO:PROJECTION:MINIMAP2_ALIGN' {
        memory = { 12.GB * task.attempt }
    }
    withName: 'FOMO:PROJECTION:MINIMAP2_INDEX' {
        memory = { 12.GB * task.attempt }
    }
}
```

Then add `-c resources_lycaenidae.config` to every command line in
`fomo_runs/lycaenidae/all_vs_all.neflow_commands.txt`, i.e.

```
nextflow run main.nf -profile crg -c /…/fomo_runs/lycaenidae/resources_lycaenidae.config --input … --outdir … -resume
```

Use the absolute path — those commands are dispatched through `send_array.sh` and the
working directory at execution time is not guaranteed to be the run directory. A `-c` file
takes precedence over profile blocks in `nextflow.config`, so this overrides the
`MINIMAP2_INDEX` `withName` in `conf/crg.config` without editing it. Changing only `memory`
leaves every task hash intact, so `-resume` still works across the shared workDir.

## 5. Deferred — bundle with the next intentional full re-run

Do **not** apply these as a hot patch; each invalidates the cache for the whole projection
stage across all cached targets in the shared workDir.

- `MINIMAP2_INDEX` cpus **8 → 2** (observed 1.15–1.54 cores). The compute saving is
  negligible (~0.03 cpu·h per target), but the scheduling saving is not: this task is a
  serialization point that everything downstream waits on, and an 8-cpu/32 GB request for
  a 14–22 s job queues far worse than a 2-cpu/12 GB one. Measured submit→done was 107 s for
  22 s of work on `P. icarus`.
- `process_low` cpus 2 → 1 for `SEQKIT_STATS` (i.e. drop the step-2 pin).

When this happens, move the `MINIMAP2_INDEX` cpus change into `conf/crg.config`'s existing
`withName` block (keep its explanatory comment, and update it — the "OOM-killed at 12 GB"
note refers to a ~3 Gb genome and must not be read as applying to Lycaenidae).

## 6. Follow-up — make minimap2 memory size-aware instead of clade-scoped

This retires step 4's per-clade file. In `conf/crg.config`, replace the static
`MINIMAP2_INDEX` memory with a closure derived from the target assembly size. Anchors from
evidence: ~0.5 Gb → 6–7 GB peak; ~3 Gb → 25–30 GB peak (per the existing comment). So
roughly `4 GB + 12 GB per Gb of assembly`.

```groovy
withName: 'FOMO:PROJECTION:MINIMAP2_(INDEX|ALIGN)' {
    memory = { (4.GB + 12.GB * Math.ceil((reference?.size() ?: 1_000_000_000) / 1e9)) * task.attempt }
}
```

Before using this, **verify the input variable name in each module** —
`modules/nf-core/minimap2/index/main.nf` and `.../align/main.nf` declare their reference
input separately and the name may differ (`reference` vs `fasta` vs a `tuple val(meta2), path(...)`
binding). A dynamic directive that references a non-existent variable fails at task
submission. Test on one Lycaenidae target **and** one large-genome target before adopting,
and keep the `?:` fallback so a missing binding degrades to 16 GB rather than crashing.

## 7. Out of scope but strongly evidence-backed: per-task scheduling overhead

Flagging because the audit surfaced it and it dwarfs everything above. Comparing `duration`
(submit → complete) against `realtime` (execution) in the traces:

| Run | Fresh (`COMPLETED`) tasks | Actual compute | Queue + poll overhead |
|---|---|---|---|
| `P. icarus` | 1424 | **1.31 h** | **16.96 h** |
| `A. agestis` | 684 | 1.38 h | 5.33 h |

(The overhead column counts `COMPLETED` rows only. `13_resource_rightsizing_audit.py`
reports the figure over *all* rows, so it shows 13.99 h for `A. agestis` — that includes the
740 `CACHED` carry-overs' original waits. Neither number is wrong; just don't compare them.)

Median per-task gap is 25–57 s against median 0.5 s of work. `conf/crg.config` sets
`executor.pollInterval = '30 sec'`, which plausibly accounts for a large share of that gap
(Nextflow detects task completion on the poll tick), with real SLURM queue time on top —
the trace cannot separate the two. Two cheap experiments, **not** part of this plan's scope:

1. Lower `executor.pollInterval` to `'5 sec'` and re-measure the gap on one target.
2. Consider whether the ~900 trivial `process_single` text-processing tasks (`*_TO_MQC`,
   `GFF_STATS*`, `FILTER_*`) should be separate SLURM jobs at all — that is a pipeline
   redesign, not a config tweak, and should be its own plan.

Raise these with the user before acting; they were not part of the resource-allocation ask.

---

## Verification protocol

Resource directives are not retroactive — a cached task keeps the trace metrics of its
original execution. So you must generate a **fresh** trace:

1. Pick one target and run it with a scratch workDir so nothing resumes:
   `nextflow run main.nf -profile crg -c …/resources_lycaenidae.config --input …/samplesheets/samplesheet_Polyommatus_icarus_265386.csv --outdir <tmp-outdir> -w <tmp-workdir>`
   Do **not** point a validation run at the shared workDir.
2. Run `python3 plans/13_resource_rightsizing_audit.py <tmp-outdir>/pipeline_info/execution_trace.txt`.

Acceptance criteria:

- [ ] `awk -F'\t' 'NR>1 && $9!="1"' execution_trace.txt` is empty — **zero retries**. Any
      row with `attempt > 1` means a tier was cut too far; check `exit` (137 = OOM,
      143 = timeout) and raise that specific label rather than reverting the whole change.
- [ ] No task has `status` other than `COMPLETED` / `CACHED`.
- [ ] Reserved memory ≤ **10 GB·h** for the run (was 21.5).
- [ ] Per-process `requested / max_peak_rss` ratio falls in **1.5×–5×** for every process
      except the `*_TO_MQC` / `GFF_STATS*` family, whose peaks (5–30 MB) are below what
      Nextflow's trace poller can measure reliably on sub-second tasks — 2 GB is already
      the practical floor there.
- [ ] The audited `gffcompare` outputs are byte-identical to the previous run for the same
      target (resource changes must not alter results — this is the real regression check).
- [ ] `-resume` on an already-completed target re-uses the cache: re-run one finished target
      and confirm the trace is ~100% `CACHED`. If it is not, a `cpus` change slipped in
      somewhere it shouldn't have (Constraint 2).

## Rollback

Steps 1–3 are three edits in two tracked files (`conf/base.config`,
`modules/local/td2_predict.nf`) — revert the commit. Step 4 is a new untracked file plus a
`-c` flag in `all_vs_all.neflow_commands.txt`; delete the file and drop the flag. Because
none of steps 1–4 change any task hash, rolling back does not invalidate the shared workDir
either.

## Appendix: measured data, for reference during implementation

Per-process request vs. peak, `P. icarus` / `A. agestis`, max across all instances. Ratios
above ~200× are on sub-second tasks where the trace poller under-samples peak RSS; the
tier is still wrong, but do not treat the exact ratio as precise.

```
process                                   n    requested        max peak RSS   max cores  max real
PROJECTION:MINIMAP2_ALIGN                 80   24 GB/4cpu/1h    6.70 / 6.80 G  3.66/3.40  89 / 82 s
PROJECTION:MINIMAP2_INDEX                  1   32 GB/8cpu/1h    6.70 / 5.90 G  1.15/1.54  22 / 14 s
PREPROCESSING:TD2_NONCODING:TD2_PREDICT   20   24 GB/4cpu/1h    0.60 / 0.60 G  1.07/1.07  15 / 15 s
PROJECTION:*:TD2_PREDICT                  20   24 GB/4cpu/1h    0.43 / 0.47 G  1.24/1.19   5 / 13 s
PREPROCESSING:EXTRACT_SEQUENCES           40   12 GB/2cpu/10m   0.18 / 0.18 G  0.95/0.95  26 / 26 s
PREPROCESSING:EXTRACT_DECOY_SEQUENCES     40   12 GB/2cpu/10m   0.15 / 0.15 G  1.00/1.00  33 / 33 s
PROJECTION:COMBINED_GTF_TO_GFF             4   12 GB/2cpu/10m   0.23 / 0.21 G  0.95/0.98   6 /  6 s
PROJECTION:EXTRACT_PROJECTED_LNC          20   12 GB/2cpu/10m   0.07 / 0.01 G  0.99/0.97  11 /  6 s
PREPROCESSING:SEQKIT_STATS                80   12 GB/2cpu/10m   0.03 / 0.03 G  0.40/0.40   2 /  2 s
PROJECTION:GFFCOMPARE_COMBINE              4    6 GB/1cpu/10m   0.76 / 0.78 G  0.78/0.72  18 / 13 s
BENCHMARKING:GFFCOMPARE                   84    6 GB/1cpu/10m   0.39 / 0.38 G  0.93/0.91   5 /  5 s
PREPROCESSING:RELOCATE_LOCI               40    6 GB/1cpu/10m   0.36 / 0.36 G  0.93/0.93   6 /  6 s
PROJECTION:SAMTOOLS_STATS                 80    6 GB/1cpu/10m   0.34 / 0.35 G  0.89/0.88   6 /  3 s
REPORTING:MULTIQC                          1    6 GB/1cpu/10m   0.30 / 0.30 G  0.14/0.25  45 / 25 s
CONSENSUS_TOP:GFFCOMPARE_TOP               2    6 GB/1cpu/10m   0.15 / 0.13 G  0.87/0.91   2 /  2 s
~880 other process_single tasks                 6 GB/1cpu/10m   <0.05 G        <0.90      <25 s
```

Caveat carried from the audit: 1216 of 1424 tasks ran under 5 s, so their `peak_rss` may be
under-sampled by Nextflow's trace poller. Every number driving a decision above was taken
from the 205–208 tasks with `realtime >= 5 s`, where measurement is reliable. Also note the
two traces are not independent samples — 740 task hashes are shared (`CACHED` carry-overs
from the same shared workDir), so `EXTRACT_SEQUENCES`, `SEQKIT_STATS`,
`PREPROCESSING:TD2_NONCODING:TD2_PREDICT` and `GUNZIP_FASTA` are literally the same
executions in both. The fresh-measurement overlap covers `MINIMAP2_*`,
`TD2_NONCODING_PROJ`, and the projection/benchmarking/consensus stages.

One more caveat for whoever touches minimap2: `MINIMAP2_ALIGN` reaches **28 GB peak vmem**
against 6.8 GB RSS (thread stacks + mmap'd index). SLURM `--mem` is RSS/cgroup-based so
12 GB is safe on `genoa64`, but do not set a virtual-memory limit below ~32 GB for this
process if a queue that enforces one is ever used.
