# Plan: CPU/memory right-sizing — no CPU escalation on retry, smaller batch and minimap2 requests

**Status: implemented** (`conf/base.config`, `conf/crg.config`). Not yet validated on a
full cluster run — see **Verification**.

Follows `plans/13` (right-sizing from Lycaenidae traces) and `plans/17` (executor tuning).
Those predate plans 15/18/21–23: `MINIMAP2_ALIGN` now aligns every source in one task per
(target, gene type), and the two `*_BATCH` processes added in plans/18 were deliberately
given the placeholder `process_high` tier "to be trimmed later, from real numbers"
(plans/18 §Decisions and step 7). This plan is that trim.

## Evidence

| Run | Where | Scale | Tasks |
|---|---|---|---|
| `curious_leavitt` | apiaceae, `--include_decoy --include_mrna` | S=5, T=35 | 1581 completed, 60.3 min |
| mollusca main | `fomo_runs/mollusca/results/pipeline_info/` | S=191, T=1 | 2501 (2470 cached) |
| mollusca OctVul | `fomo_runs/mollusca/results_OctVul/pipeline_info/` | S=15, T=1 | 627 (2 failed, retried) |

Caveats:
- The mollusca main trace **predates plans/20** (it still has `FILTER_LNC_GFF`) and ran with
  the one-off `mollusca/index_highmem.config` (INDEX 256 GB, ALIGN 128 GB) for the 14.7 Gb
  *Vampyroteuthis* target. Its memory figures are not pipeline defaults, and its CPU figures
  describe one giant genome.
- "Utilisation" = trace `%cpu` / requested cpus. `%cpu` is a whole-task average, so it hides
  phases (e.g. multi-threaded minimap2 followed by a single-threaded `samtools sort`).

### Measured use vs request

| Process | Request (before) | Peak RSS (max) | Cores used (avg) |
|---|---|---|---|
| `MINIMAP2_ALIGN` | 4 cpu / 24 GB | 20.3 GB apiaceae, 15.4 GB OctVul | ~1.0 apiaceae, 0.8 OctVul, 2.65 on the 14.7 Gb genome |
| `MINIMAP2_INDEX` | 8 cpu / 32 GB | 28.3 GB apiaceae | ~1.0 apiaceae, 0.8 OctVul, 2.1 on the 14.7 Gb genome |
| `GFF_STATS_PROJECTED_BATCH` | 12 cpu / 72 GB | 0.15 GB | ~1.1 |
| `GFFCOMPARE_BATCH` | 12 cpu / 72 GB | 0.53–0.73 GB | 1.2–2.9 |
| `SEQKIT_STATS` | 2 cpu / 4 GB | 0.03 GB | ~0.1–0.2 |

Apiaceae `MINIMAP2_ALIGN` alone is 37.6 of ~44 requested CPU-h, ~75 % of it idle.

### Queue wait (apiaceae timeline, submit → start)

| Process | median | max |
|---|---|---|
| `MINIMAP2_ALIGN` | 892 s | 1759 s |
| `GFF_STATS_PROJECTED_BATCH` / `GFFCOMPARE_BATCH` | ~438 s | 745 s |
| `MINIMAP2_INDEX` | 331 s | 607 s |
| 1-cpu tasks | ~5 s | — |

The long waits are on the large requests, so smaller requests should also shorten wall
time. That is an expectation, not a measurement — SLURM QOS/fairshare was not inspected.

### Retries ask for more CPUs

`conf/base.config` said "`cpus` deliberately does NOT escalate", but three tiers did:
the process default `{ 1 * task.attempt }`, `process_medium` `{ 4 * task.attempt }` and
`process_high` `{ 12 * task.attempt }`. An OOM retry of `MINIMAP2_ALIGN` therefore asked for
8, then 12 cpus, and a batch retry for 24/36 cpus with 144/216 GB. Extra cores never fix an
OOM or a time-limit kill, and a bigger request only waits longer in the queue.

## Changes

1. **CPUs never escalate on retry.** The process default, `process_medium` and
   `process_high` get fixed `cpus` (1, 4, 12). Memory and time keep their `* task.attempt`
   escalation, which is what actually recovers OOM and time-limit kills.
2. **The two batch processes** (`GFF_STATS_PROJECTED_BATCH`, `GFFCOMPARE_BATCH`): a
   `withName` block in `conf/base.config` sets `cpus = 4`, `memory = { 4.GB * task.attempt }`.
   Peak is ≤ 0.73 GB with 12 workers; 4 workers keep the `xargs -P` loop parallel for large
   S. `time` stays at `process_high`'s 2 h. `process_high` itself is left as a generic tier
   (its "unused" comment was stale and is corrected).
3. **minimap2 CPUs.**
   - `MINIMAP2_ALIGN`: `cpus = 2` (was 4), in the existing `conf/base.config` block.
     `-t 2` lets minimap2 run up to 3 threads (it adds one I/O thread) and gives
     `samtools sort` `-@ 1`. Memory is unchanged (85 % of the request on apiaceae).
   - `MINIMAP2_INDEX`: `cpus = 3` (was 8), in `conf/crg.config`. minimap2 uses **at most three
     threads when indexing**, whatever `-t` says (minimap2 man page, `-t`), which matches the
     2.1 cores measured on the 14.7 Gb genome. Memory unchanged.
4. **`SEQKIT_STATS`**: drop its `cpus = 2` override and fall back to `process_low`'s 1 cpu.
   plans/13 deferred this "until the next intentional full re-run", which change 3 forces.

### Not changed, and why

- **`process_low` memory stays at 4 GB.** It also covers `TD2_PREDICT`, which peaked at
  1.4–14.2 GB in the mollusca trace (42 of 192 tasks above 1 GB). Lowering it would buy
  little — memory for 1-cpu jobs rarely limits scheduling.
- **`process_single` (1 cpu / 2 GB)** is already sized.
- **Local executor tier** — dropped. Short tasks wait ~5 s in the queue; the rest of their
  ~33 s duration is container start, NFS staging and exit detection, which a local task
  pays too. Too little gain for the head-node risk.
- **Dropping `MINIMAP2_INDEX` when only lncRNA runs** — dropped. Two indexing paths
  depending on `--include_mrna`/`--include_decoy` complicate the code for one job per target.
- **`EXTRACT_SEQUENCES` time** — unchanged. The OctVul failures (killed at 14m 36s, then
  2m 55s on retry) look like a stalled first attempt, which a longer limit would only prolong.

## Consequences for `-resume`

`task.cpus` is interpolated into the command of `MINIMAP2_INDEX` (`-t`), `MINIMAP2_ALIGN`
(`-t`, `samtools sort -@`), `SEQKIT_STATS` (`--threads`) and both batch processes
(`xargs -P`), so changes 2–4 change their task hashes. Resuming an existing run re-runs
every index and alignment, and therefore everything downstream of projection, plus SeqKit
stats and their MultiQC adapters. Memory-only changes do not re-hash. Change 1 does not re-hash
anything on its own (attempt 1 resolves to the same value).

## Verification

- Done: `nextflow config -profile crg` parses, and the resolved process block carries the
  new values.
- To do: next full run (e.g. apiaceae, fresh or resumed). From `execution_trace.txt`:
  - requested `cpus`/`memory` match the table above, including on retried rows (no
    escalating cpus);
  - no new exit 137/140/143 rows on the batch processes;
  - `MINIMAP2_ALIGN` / `MINIMAP2_INDEX` `realtime` vs `curious_leavitt` (alignment may be
    slower per task — plans/13 measured 3.4–3.66 cores on small Lycaenidae genomes before
    plans/15);
  - queue wait (start − submit) for ALIGN, INDEX and the batch processes.

## Open follow-ups (not in this change)

- `TD2_PREDICT` peaked at 14.2 GB on a 4 GB request without being killed. Check whether
  SLURM enforces memory on genoa64, and whether TD2 needs its own tier.
- Profile one `EXTRACT_SEQUENCES` task (3 % CPU at S=191; likely NFS I/O on the
  uncompressed genome).
- Fusing the short per-target chain (`FILTER_ALLMODELS → BAM_TO_GFF → SPLIT_GFF_BY_SOURCE`,
  `SAMTOOLS_STATS → SAMTOOLS_TO_MQC`) is the remaining lever on per-task overhead; it touches
  the `ext.prefix`/`publishDir` wiring and needs its own plan.
