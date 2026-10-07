# 19 — RUN_SUMMARY_TABLES at scale: why it timed out, and what to do about it

## Problem, with numbers

`FOMO:RUN_SUMMARY:RUN_SUMMARY_TABLES` (`bin/run_summary_tables.py`) is the pipeline's one
remaining serial bottleneck: a single task, no target dimension, that digests the whole
run's `mqc_files` union into the run-level report. On the Lepidoptera run it was bumped to
`2h × task.attempt` / `8 GB × task.attempt` (`conf/crg.config`, following the earlier
`RUN_SUMMARY_TABLES (run_summary)` walltime incident) — and it **still timed out**:

```
28387428   TIMEOUT   02:00:26   ReqMem=8G   MaxRSS=3.5GB
```

**Not memory** — `MaxRSS` was under half the request. Pure throughput. The task's own work
dir, now that the attempt has run to completion of its inputs, shows exactly why:

```
613,991 files in mqc/        (427,720 *_gffstats_projection_mqc.tsv, 180,512 *.gffcompare.stats,
                               ~5,300 other *_mqc.tsv, ~1,260 SAMTOOLS_STATS *.raw.stats)
531     files in top_sources/
```

Three serial, single-threaded passes touch nearly all of them:

1. **`Inputs.__init__`, `bin/run_summary_tables.py:138`**:
   `files = sorted(p for p in mqc_dir.rglob("*") if p.is_file())` — `is_file()` **follows
   the symlink** (everything here is a Nextflow-staged symlink) and `stat()`s the target
   over NFS. **~614,000 stat round-trips**, before a single byte of content is read.
2. **`read_tsvs()` via `pick()`**: opens and `csv.DictReader`s every `*_mqc.tsv` it picked
   out of that list — **~432,700 file opens**.
3. **`Inputs.accuracy()`**: opens and parses every `.gffcompare.stats` file individually via
   `fomo_stats.parse_accuracy` — **~180,500 more file opens**.

**≈1.23 million individual file-touching syscalls, serially, in one Python process, over
NFS.** At even a modest 5–10 ms average round-trip under this run's concurrent cluster load,
that alone is 1.7–3.4 hours — before any parsing or table-building. This is not a script bug
in the sense of wrong logic (every function here is correct and was verified end-to-end on
`-profile test`); it simply was never exercised at real full scale until now. The plans/18
batching cut *task count* (the OOM problem); it deliberately did **not** cut *file count*,
since per-source granularity is what the run-level report is for.

## Option 1 — Optimise the script

Two independent, additive changes, both low-risk (pure engineering, no output change):

### 1a. Delete the stat call in the listing (free, zero risk)

`mqc/` is populated **entirely** by Nextflow-staged symlinks with no subdirectories — every
entry is, by construction, a file. `is_file()` is the wrong tool twice over: it follows the
symlink (an NFS round-trip to a *different* directory) when `entry.is_symlink()` alone
(a local `lstat` on the entry itself, no round-trip to the target) already answers "is this
a real entry to keep". Replace:

```python
files = sorted(p for p in mqc_dir.rglob("*") if p.is_file())
```

with an `os.scandir()`-based version that only asks about the entry itself:

```python
with os.scandir(mqc_dir) as it:
    files = sorted(Path(e.path) for e in it if e.is_symlink() or e.is_file(follow_symlinks=False))
```

Removes **~614,000 stat round-trips** for zero behavioural change. Keep
`_assert_unique` as-is — it's an in-memory dict pass, negligible cost.

### 1b. Parallelise the two per-file read passes (the big win)

`read_tsvs()` and `Inputs.accuracy()` are two independent loops that each just
open-read-parse one small file at a time — the textbook case for
`concurrent.futures.ThreadPoolExecutor`: Python releases the GIL during blocking I/O, so a
thread pool gets real wall-clock parallelism here with no multiprocessing complexity, exactly
the same "many independent small I/O ops" shape `xargs -P` already solved for
`GFF_STATS_PROJECTED_BATCH`/`GFFCOMPARE_BATCH` (plans/18) — just inside Python instead of
bash this time, because these are read-only single-process operations with no external
binary to invoke per file.

```python
from concurrent.futures import ThreadPoolExecutor

def read_tsvs(paths: Sequence[Path], workers: int) -> List[dict]:
    def read_one(path: Path) -> List[dict]:
        with path.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for row in rows:
            row["__file"] = path.name
        return rows

    rows: List[dict] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for chunk in pool.map(read_one, paths):
            rows.extend(chunk)
    return rows
```

Same shape for `Inputs.accuracy()`'s loop over `self.stats_files`. Add an `--io-workers`
CLI flag (default e.g. 32, wired from `task.cpus` via `conf/modules.config`'s
`RUN_SUMMARY_TABLES` block so it scales with whatever the process is actually granted —
`process_single`'s 1 cpu is the wrong ceiling for a thread pool doing I/O-bound work, so this
process's `cpus` should also go up, not just its `memory`/`time`). Order of `rows`/`acc`
matters nowhere downstream (both are immediately re-indexed into dicts keyed by identity —
`annotation_index()`, `accuracy()`'s own return, etc. — never consumed positionally), so
`pool.map`'s result order (which happens to be preserved anyway) needs no special handling.

**Combined effect**: 1a removes ~614K round-trips outright; 1b turns the remaining ~613K
serial opens into ~613K/`workers` wall-clock opens. At 32 workers that is nominally a ~20–30×
speedup on the I/O-bound portion — plausibly turning "didn't finish in 2h, projected 3–4h+"
into low tens of minutes. **Verify this projection on the real numbers** (Verification
section) before trusting it for the next full run; NFS may not scale linearly past some
concurrency, so the flag being tunable matters.

**Risk**: low. No logic changes, no output-shape changes, fully covered by diffing against
a known-good run (`-profile test`, and — uniquely available right now — the *actual* stuck
run's already-fully-materialized `mqc/`/`top_sources/` directories, see Verification).

## Option 2 — Run the remaining step(s) outside Nextflow (analysis, as requested)

### Why it's actually easy right now, specifically

`RUN_SUMMARY_TABLES` and `MULTIQC_RUN_SUMMARY` are the **last two nodes in the entire DAG** —
nothing inside Nextflow consumes their output. Everything upstream for this run is already
done: `WorkflowStats` from the failed session reports `succeededCount=1251`,
`cachedCount=33580`, and every per-target `results/targets/<t>/` tree is already published.
The failed task's own work dir already holds:

- a **fully-formed `.command.sh`** (Nextflow already resolved the exact file list and CLI
  flags — nothing to reconstruct)
- `mqc/` and `top_sources/` **fully staged and complete** (613,991 + 531 symlinks — this is
  not partial data left over from a crash, it's the complete input set)

So concretely: `sbatch` that same `.command.sh` (or a copy with a generous, one-off `-t`,
e.g. 8h) directly against that work dir, let it run unsupervised to completion, then hand-run
the equivalent `MULTIQC_RUN_SUMMARY` invocation (also trivial to reconstruct — it's a single
MultiQC call over the just-produced tables plus the two config YAMLs), and copy the results
into `results/summary/tables/` and `results/summary/multiqc/` — the same paths
`conf/modules.config`'s publishDir rules would have used.

### Trade-offs

**Pros:**
- Unblocks *this* run today, independent of any code fix or further Nextflow attempts.
- Removes the fragile, repeatedly-OOMing Nextflow head process from the critical path for
  the slowest remaining step — the head only needs to poll a `squeue`/file-existence check,
  not hold the DAG open for hours.
- Zero new code to trust before getting a result.

**Cons / risks:**
- **Nextflow's own bookkeeping goes stale.** `-resume`'s cache DB, `execution_trace.txt`,
  the timeline/report HTML will all show these two tasks as failed/missing even though the
  outputs exist on disk. A later `-resume` of this run UUID will try to re-run them — harmless
  (cheap no-op if the script's fixed and fast, a repeat of the same problem if it isn't) but
  worth knowing about, not hand-patching the cache DB is the simpler choice here since this
  is the terminal step.
- **Manual publish must exactly match** the paths/filenames `conf/modules.config`'s
  `'.*:RUN_SUMMARY_TABLES$'` / `'.*:MULTIQC_RUN_SUMMARY$'` rules use, or anyone (person or
  tool) expecting the "official" `results/summary/` location won't find it.
- **One-off, not a fix.** Doesn't help a future full rerun, a bigger clade, or anyone else
  hitting this same step — it only rescues *this* run's already-computed upstream work.
- In tension with reproducibility (a hand-run step outside the pipeline's own provenance
  trail) — acceptable as a one-time rescue of 90%-complete production work, not as standing
  practice.

**Verdict:** good as an **immediate, parallel unblock** for the current run — it costs
nothing to prepare and doesn't compete with fixing the script. Not a substitute for Option 1
for anything that runs again.

## Option 3 — Other alternatives considered

### 3a. Structural: cut the file *count*, not just the read speed

Have `GFF_STATS_PROJECTED_BATCH`/`GFFCOMPARE_BATCH` (plans/18) each also emit one
concatenated per-target digest (rolling up their S-many internal per-source rows into one
file) instead of routing every individual per-source file into the run-wide `mqc_files`
union. Would cut `RUN_SUMMARY_TABLES`'s input count by roughly a factor of S (~500×): from
~614K down to ~T-scale (~1,200). The most durable fix in principle — same "batch small files
into fewer" philosophy as plans/18 — but it's real surgery on modules already hardened once
(the naming-grammar and MultiQC-routing correctness plans/18 needed care for), and it trades
away the ability for MultiQC/ad-hoc inspection to reach an individual source's raw row
without first computing the digest. **Worth doing only if 1a+1b don't scale far enough** for
a meaningfully bigger future run — no evidence yet that they won't.

### 3b. Nextflow-native map/reduce

Split the `mqc_files` union into N chunks (channel-level `.collate()`/round-robin by
species), run N small "partial digest" tasks in parallel as ordinary SLURM jobs, merge with
one final light task. Gets the same parallel-I/O win as 1b but spread across the cluster
instead of one long-lived task, and doesn't need one big/long single-task SLURM allocation
babysat by the fragile head process for hours. More invasive than 1b (new subworkflow, a
real merge step) and probably unnecessary if in-process threading (1b) gets this under
~20–30 minutes — keep in reserve if NFS turns out to cap per-client concurrency lower than
expected.

### 3c. Stop-gap: just raise `time` further

Zero engineering risk, but doesn't fix anything, and the ~1.2M-serial-touch estimate above
means even 6–8h isn't a *safe* margin, only a *probably-enough* one — and it keeps the
already-OOM-prone head process alive and babysitting for that whole window. **Not
recommended alone.**

## Recommendation — sequenced

1. **Now, to unblock this run today:** Option 2. Everything needed already exists on disk;
   this costs nothing and runs in parallel with step 2. I have not run anything yet — this
   needs your go-ahead since it's a direct action against the production run/work dir.
2. **Before any future full rerun:** implement Option 1 (1a + 1b) in
   `bin/run_summary_tables.py`, plus the matching `cpus`/`ext.args` wiring in
   `conf/modules.config`/`conf/crg.config`. Verify on `-profile test`, **and** — a uniquely
   good, free, full-scale integration test available right now — point the fixed script
   directly at this stuck run's already-materialized `mqc/`/`top_sources/` directories and
   diff its output against whatever Option 2's manual run produces. That proves correctness
   *and* gives a real measured speedup number before trusting it for the next run.
3. **Only if a future run is meaningfully bigger and 1a+1b still isn't comfortably under a
   couple hours:** revisit 3a or 3b.

## Verification (for step 2, once implemented)

1. `-profile test,singularity` end-to-end run, diff `run_summary.json` and every
   `run_*_mqc.{tsv,json}` byte-for-byte against current `main` — this is a pure performance
   change, zero output diff expected.
2. Real-scale check against the stuck run's own `mqc/`/`top_sources/` (read-only — copy
   the work dir or point `--mqc-dir`/`--top-sources-dir` at it, don't run in place): measure
   wall time at a few `--io-workers` values, confirm it lands comfortably inside the
   `conf/crg.config` `time` budget, and diff its `run_summary.json` against Option 2's
   manually-produced one — should be identical (same inputs, same logic, only the read
   strategy changed).
