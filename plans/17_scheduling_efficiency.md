# Plan: Scheduling efficiency — executor tuning + removing redundant gunzip tasks

## Context

The Lepidoptera run `agitated_legentil` (Aug 20 2026, 15:58:15 → 16:59:53) aborted after
61.6 min with `java.io.IOException: Disk quota exceeded` while writing a task wrapper into
`workDir = /nfs/scratch01/rg/fzanarello/work_fomo` (`.nextflow.log:14149`). The samplesheet
has 1260 rows (531 `both`, 729 `target`) → S=531, T=1260, T_g=531. It never left
preprocessing: 200 of 1260 `MINIMAP2_INDEX` done, zero `MINIMAP2_ALIGN`.

Auditing `results/pipeline_info/execution_trace.txt` (4119 `COMPLETED`, zero failures):

| | median | mean | total |
|---|---|---|---|
| `duration` (Nextflow-observed) | 79.0 s | 77.1 s | **88.2 h** |
| `realtime` (the script itself) | 4.8 s | 14.9 s | **17.0 h** |
| overhead | 60.0 s | 62.3 s | **71.2 h** |

**Only 19.3% of allocated task time did useful work.** The median task did 4.8 s of work
inside a 79 s slot, and 66% of all tasks (2738) ran under 15 seconds. Splitting the 60 s
of overhead from the log (submit → SLURM exit → Nextflow notices): **~25 s is Nextflow's
poll loop**, ~36 s is SLURM queue wait + singularity/stage-in.

Two facts establish that Nextflow, not SLURM, was the binding constraint:

- `.nextflow.log:75` — `Creating task monitor for executor 'slurm' > capacity: 100`.
  `queueSize` is unset anywhere in the repo, so it is the SLURM default of 100. In-flight
  jobs (Nextflow's view) sat at median 77, max 100.
- `sacct` shows the run's jobs auto-assigned to QOS **`test`** (10-min `-t`) and **`vshort`**
  (1 h `-t`), each with `MaxJobsPU = 96` and no CPU/memory `MaxTRESPU` cap. The association
  sets no `MaxJobs`/`MaxSubmit`. So SLURM would have run ~96 concurrent jobs *per QOS tier*
  while Nextflow was holding a 100-slot window with a third of each slot's life spent
  waiting for a 30 s poll.

Submission itself is not the bottleneck: the median inter-submission gap is 0.196 s
(≈5 tasks/s of sbatch capacity), while the run averaged 1.14 tasks/s. The gap was
slot starvation, not sbatch cost.

The overhead is concentrated in trivially short tasks, which is what makes this fixable:

| process | n | med realtime | med duration | Σ overhead | med %cpu |
|---|---|---|---|---|---|
| `PROJECTION:GUNZIP_TARGET` | 1151 | 21.5 s | 87.0 s | **20.68 h** | 15.3% |
| `BENCHMARKING:GUNZIP_TARGET_GFF` | 445 | **1.1 s** | 73.0 s | **8.04 h** | 22.3% |
| `PREPROCESSING:GUNZIP_FASTA` | 411 | 15.7 s | 86.0 s | **7.79 h** | 19.6% |

This plan addresses the two cheap levers: **(1)** executor tuning, **(2)** deleting the
gunzip tasks whose consumers provably read `.gz` natively. Both reduce the per-task latency
that dominates the run; (2) also reduces what gets written into the over-quota work dir.

**Explicitly out of scope** (needs its own plan): the quadratic task counts that make a
1260-species run infeasible regardless of latency — `GFF_STATS_PROJECTED` at S×T = 669,060
tasks and the benchmark `GFFCOMPARE` at `(F·S·D + 2·F·D)·T_g` = 283,023. Those are 952k of
the ~968k tasks a full run needs. This plan does not make a 1260-species run finish; it
makes each task ~3× cheaper and unblocks the disk quota.

---

## Part 1 — Executor tuning (`conf/crg.config` only)

Scoped to the `crg` profile so `-profile test` and local dev are untouched. The profile
already has the only executor-scope setting in the repo (`conf/crg.config:60`,
`executor.pollInterval = '30 sec'`); this replaces that line with a full `executor` block.

### 1.1 Replace the bare `executor.pollInterval` line

In `conf/crg.config`, inside the `crg` profile, next to `singularity.enabled` /
`singularity.autoMounts`, replace:

```groovy
        executor.pollInterval = '30 sec'
```

with:

```groovy
        // ── Executor tuning ──────────────────────────────────────────────────
        // Measured from run `agitated_legentil` (see plans/17): the median task
        // did 4.8 s of work inside a 79 s slot. ~25 s of that was this poll loop
        // and the rest SLURM queue wait — 80.7% of allocated task time was
        // overhead. 66% of tasks run under 15 s, so a 30 s poll is sized for a
        // workload this pipeline does not have.
        executor {
            // A finished 1-second job used to sit undetected for a mean of 31 s,
            // holding one of the queueSize slots the whole time.
            pollInterval = '5 sec'

            // Default is 100 and it was the binding constraint: in-flight sat at
            // median 77 / max 100 while QOS `test` and `vshort` (which the 10-min
            // and 1-h `-t` values auto-select) each allow MaxJobsPU = 96 with no
            // CPU/mem MaxTRESPU cap, and the association sets no MaxSubmit. A
            // backlog larger than the running cap is intentional — it is what
            // refills a slot the instant one frees.
            //
            // Deliberately NOT larger: RUNNING jobs are capped SLURM-side at
            // MaxJobsPU = 96 per QOS tier (~200-300 across the tiers this pipeline
            // touches), so a deeper queue only adds PENDING jobs, sbatch ramp time
            // (~0.2 s each, single-threaded) and queue depth for genoa64's
            // scheduler to walk (default_queue_depth=500, bf_max_job_user=20).
            // 500 is ~2x the running cap: enough backlog to refill a freed slot
            // instantly, without flooding the queue.
            queueSize = 500

            // sbatch costs ~0.2 s from the single submitter thread (~5 tasks/s of
            // headroom), so no rate limit is needed; left unset deliberately.

            dumpInterval = '1 min'
        }
```

### 1.2 Guard the one fat process against the new concurrency

`conf/crg.config:20-24` gives `MINIMAP2_INDEX` 8 cpus / 32 GB. At the old 100-slot window
that was self-limiting; at 500 it is not. Partition `genoa64` has `TotalCPUs=896`
(`TRES=cpu=868,mem=4953396M`) across 7 nodes. Add a `maxForks` to that existing
`withName` block, set to **96** — a deliberate throughput-first choice: 96 concurrent
indexes claim 768 of 868 cpus and 3.0 of 4.95 TB, i.e. close to the whole partition. It
also coincides with the SLURM-side ceiling, since the 1 h `time` auto-selects QOS `vshort`
(`MaxJobsPU = 96`), so the cap mostly stops Nextflow queueing work SLURM will not run.
`MINIMAP2_ALIGN` (`time = 1.h` in `base.config`) lands in that same QOS and competes for
the same 96 slots. At T = 1260 indexes of ~68 s each this is ~15 min of indexing versus
~60 min at a more conservative 24; lower it if partition citizenship matters more than
wall clock.

```groovy
            withName: 'FOMO:PROJECTION:MINIMAP2_INDEX' {
                cpus   = { 8     * task.attempt }
                memory = { 32.GB * task.attempt }
                time   = { 1.h   * task.attempt }
                // 8 cpu / 32 GB each, vs process_single (1 cpu/2 GB) elsewhere.
                // 96 matches QOS vshort's MaxJobsPU, i.e. as fast as SLURM will
                // run them; drop to ~24 to leave the partition room.
                maxForks = 96
            }
```

**Do not add `maxForks` in `conf/base.config`.** Per the trap documented at
`conf/base.config:80-92`, `crg.config` is included *before* `base.config`, and when two
`withName` selectors match, the later declaration wins regardless of specificity — a broad
selector in `base.config` would silently override this one.

### 1.3 Leave alone

- `resourceLimits` (`conf/crg.config:12-16`) — unchanged. The auto-selected QOS follows the
  requested `-t`, so shorter walltimes already land in the higher-`MaxJobsPU` tiers by
  themselves. Do **not** pin a QOS via `clusterOptions`: `normal` caps at
  `MaxJobsPU = 32`, three times worse than what the run already got.
- `exitReadTimeout` — default 270 s is right for NFS with `cache = 'lenient'`.

---

## Part 2 — Delete the gunzip tasks whose consumers read `.gz` natively

`modules/local/maybe_gunzip.nf` (created by `plans/08_mixed_compression_inputs.md`) is
aliased three times. Two of its premises have since gone stale: that plan's line 30-32 says
*"gffread, minimap2 (index/align), samtools faidx all need an uncompressed FASTA; AGAT needs
an uncompressed GFF3"* — but **minimap2 reads gzip natively** and **AGAT is gone**, replaced
by `gff-feature-stats`, which reads `.gff3.gz` directly (`modules/local/gff_stats.nf:5-8`).
That same plan already lost a fourth alias, `GUNZIP_RAW_SOURCE_GFF`, for exactly this reason
— the raw source GFF3 now goes into `GFF_STATS` still gzipped (`preprocessing.nf:202`).
That removal is the precedent for 2.1.

Verified consumer-by-consumer:

| Alias | Consumer | Reads `.gz`? |
|---|---|---|
| `GUNZIP_TARGET_GFF` | `FILTER_TARGET` | **Yes** — inline `if [[ "$gff3" == *.gz ]]` at `modules/local/filter_transcript.nf:20-24` |
| `GUNZIP_TARGET_GFF` | `GFF_STATS_TARGET` | **Yes** — native, `modules/local/gff_stats.nf:53` |
| `GUNZIP_TARGET` | `MINIMAP2_INDEX` | **Yes** — minimap2 gzopen's its input |
| `GUNZIP_TARGET` | `ch_ref` FASTA slot | **Dead value** — discarded as `_fa` at `projection.nf:76` |
| `GUNZIP_TARGET` | `EXTRACT_PROJECTED_LNC` (gffread `-g`) | **No** |
| `GUNZIP_FASTA` | `EXTRACT_SEQUENCES` / `EXTRACT_DECOY_SEQUENCES` (gffread `-g`) | **No** |
| `GUNZIP_FASTA` | `SAMTOOLS_FAIDX` | bgzip only — see the trap below |

### 2.1 Delete `GUNZIP_TARGET_GFF` entirely — `subworkflows/local/benchmarking.nf`

Both consumers are proven `.gz`-capable, and both derive their output names from
`ext.prefix` (`conf/modules.config:109-116` builds `<id>.target.raw` from `meta`, and
`FILTER_TARGET`'s prefix is `${meta.id}.${meta.feature_type}`), so **no published filename
changes**. This removes T_g tasks — 445 observed, 531 at full scale, 8.04 h of overhead.

1. Delete the include at `benchmarking.nf:1`.
2. Delete the `GUNZIP_TARGET_GFF(...)` call at `benchmarking.nf:27-29`.
3. Repoint the two consumers at the raw samplesheet GFF3 already carried by `ch_target_gff`
   (`benchmarking.nf:25`). `ch_target_gff` is `[meta, fasta, gff3]`, so map it to
   `[meta, gff3]` once and reuse:
   - `benchmarking.nf:32` — `ch_target_filter_in` currently `flatMap`s over
     `GUNZIP_TARGET_GFF.out.gunzip`.
   - `benchmarking.nf:71` — `ch_target_stats_gff` currently starts from the same channel.

   Read the channel twice (Nextflow forks it) or bind it to a local, matching the existing
   style. Note `FILTER_TARGET` and `GFF_STATS_TARGET` will now each decompress
   independently, so the ~1.1 s of gunzip work happens twice instead of once — a trade the
   numbers overwhelmingly justify (1.1 s of duplicated work against 73 s of scheduling
   latency plus a full-file `cp -L`).
4. Delete the now-dead selector at `conf/modules.config:295-297`. A stale `withName` is
   silently ignored by Nextflow rather than erroring, so it would otherwise just rot.
5. Update the stale doc comment at `benchmarking.nf:9` (`gff3[.gz]` wording is already
   right; check the surrounding prose still describes the flow).

**Precondition already satisfied:** `gff-feature-stats` requires each seqid's records to be
contiguous above its 200k-line batch threshold (`modules/local/gff_stats.nf:18-23`). The
raw Ensembl GFF3 is grouped, which `preprocessing.nf:202` already relies on for the
identical file class.

### 2.2 Keep `GUNZIP_TARGET` and `GUNZIP_FASTA` — and record why

gffread 0.12.7 (`modules/nf-core/gffread/main.nf`) opens the `-g` genome through a
random-access FASTA reader with no gzip or bgzf path, so `EXTRACT_SEQUENCES`,
`EXTRACT_DECOY_SEQUENCES` and `EXTRACT_PROJECTED_LNC` genuinely require plain FASTA. Both
aliases stay.

**The `SAMTOOLS_FAIDX` shortcut is a trap — do not take it.** `samtools faidx` accepts a
*bgzip* FASTA but hard-errors on plain gzip. The committed test data is bgzip
(`bin/subsample_test_data.sh:65-76` pipes the FASTA through `bgzip` and the GFF3 through
`gzip`), but the real inputs are **plain gzip** — verified on the Lepidoptera samplesheet:
`GCA_964187985.1_ilAbrSylv1.hap1.1_genomic.fna.gz` has gzip `FLG=0x00`, no `BC` extra
subfield. `assets/schema_input.json:22-35` only pattern-matches the `.gz` suffix, so nothing
enforces BGZF. A faidx shortcut would pass `-profile test` and fail on the cluster.

Add a short comment at `preprocessing.nf:30` and `projection.nf:25` naming gffread as the
sole reason each alias survives, so this is not re-litigated. Note the disk consequence
honestly: ~T decompressed genomes still land in the work dir, so **Part 2 does not by itself
solve the quota** — it removes T_g of the ~1600 gunzip tasks and their `cp -L`, and the
quota headroom comes from the user's cleanup plus the smaller task count.

### 2.3 Drop the dead FASTA join — `subworkflows/local/projection.nf:36-38`

```groovy
    ch_ref = MINIMAP2_INDEX.out.index
        .join(GUNZIP_TARGET.out.gunzip)
        .map { meta, mmi, fa -> tuple(meta.id, mmi, fa) }
```

The FASTA is discarded at `projection.nf:76` (`{ meta, reads, tid, mmi, _fa -> }`) —
`MINIMAP2_ALIGN` receives only the `.mmi`. Reduce to:

```groovy
    ch_ref = MINIMAP2_INDEX.out.index.map { meta, mmi -> tuple(meta.id, mmi) }
```

and drop `_fa` from the `multiMap` closure at `projection.nf:76`. Pure simplification: it
also removes a false ordering dependency that made every alignment wait on the target
genome being decompressed. Update the `projection.nf:33-35` comment, which currently
explains a join that no longer exists.

### 2.4 Feed `MINIMAP2_INDEX` the raw `.gz` — `subworkflows/local/projection.nf:31`

minimap2 reads gzipped FASTA, so change

```groovy
    MINIMAP2_INDEX(GUNZIP_TARGET.out.gunzip)
```

to take the samplesheet FASTA directly from `ch_target`. With 2.3 done, `GUNZIP_TARGET`
then feeds only `ch_target_fa` → `EXTRACT_PROJECTED_LNC`, which sits downstream of the
entire alignment chain — so indexing no longer waits on decompression, and gunzip tasks stop
competing for slots during the align phase.

**No task-count reduction here** (`EXTRACT_PROJECTED_LNC` still needs one plain genome per
target); the win is critical-path only. Two side effects, both benign:
- The index filename changes. `modules/nf-core/minimap2/index/main.nf:23-27` names its output
  `${fasta.baseName}.mmi` and **never reads `task.ext.prefix`**, so
  `conf/modules.config:188`'s `ext.prefix = { "${meta.id}.splice" }` is already dead config
  today; `<id>.mmi` becomes `<GCA_accession>.fna.mmi`. Nothing keys off the `.mmi` name —
  `MINIMAP2_ALIGN` takes it as a path. The `ext.args = '-x splice'` on the same block is
  load-bearing and must stay.
- The `-x splice` index of a gzipped genome is byte-identical to that of the decompressed
  one; only decompression location moves.

If this step proves troublesome, it reverts independently of 2.1–2.3.

---

## Files touched

| File | Change |
|---|---|
| `conf/crg.config` | replace `executor.pollInterval` line with an `executor {}` block; add `maxForks` to the existing `MINIMAP2_INDEX` `withName` |
| `subworkflows/local/benchmarking.nf` | delete `GUNZIP_TARGET_GFF` include + call; repoint 2 consumers at the raw gff3 |
| `subworkflows/local/projection.nf` | simplify `ch_ref`; drop `_fa`; repoint `MINIMAP2_INDEX`; comments |
| `subworkflows/local/preprocessing.nf` | comment only (why `GUNZIP_FASTA` stays) |
| `conf/modules.config` | delete the dead `GUNZIP_TARGET_GFF` selector (`:295-297`) |
| `CLAUDE.md` | task-count table (`F·D·T` gunzip line), and the `gff-feature-stats` reads-`.gz` note now applies to the target side too |
| `plans/08_mixed_compression_inputs.md` | note that the minimap2 and AGAT premises are superseded by plans/17 |

`modules/local/maybe_gunzip.nf` itself is **not** modified — two aliases still use it.

## Verification

There is no nf-test harness, no `tests/`, and no CI in this repo, so the profile smoke test
is the whole regression net — run it before and after and diff.

1. **Baseline, then after:** `nextflow run . -profile test,docker` on the committed
   `assets/samplesheet.csv` (5 species, one contig each). Must complete. Diff the two
   `results/` trees — **published filenames and GFF/stats content must be byte-identical**
   except for the `.mmi` rename from 2.4. This is the key assertion: every affected
   `ext.prefix` is meta-derived, so nothing user-visible should move.
2. **Confirm the tasks are gone:** `grep -c GUNZIP_TARGET_GFF results/pipeline_info/execution_trace.txt`
   → 0, and the `GUNZIP_TARGET` row count is unchanged (it still serves gffread).
3. **Confirm the executor block took:** the new run's `.nextflow.log` must show
   `capacity: 500; pollInterval: 5s` in the `Creating task monitor` line — that log line is
   the only direct proof the settings were parsed.
4. **Plain-gzip guard — partly covered by the test profile, contrary to first
   assumption.** `bin/subsample_test_data.sh` pipes the FASTA through `bgzip` but the GFF3
   through plain `gzip`, and this was verified on the committed files: all five
   `*.gff3.gz` are **plain gzip**, all five `*.fna.gz` are **BGZF**. So `-profile test`
   *does* exercise 2.1 (plain-gzip GFF3 straight into `FILTER_TARGET` and
   `GFF_STATS_TARGET`) — the consequential half. What it does **not** cover is 2.4:
   minimap2 indexing a *plain-gzip* assembly, since the test FASTAs are BGZF. BGZF is
   gzip-compatible and minimap2 gzopen's its input, so this is low risk, but to close it
   run 3–4 rows of the real Lepidoptera sheet (at least one `both` and one `target`) under
   `-profile crg`. That real-data run is also the only thing that proves no faidx shortcut
   crept back in.
5. **Measure the win** on that same small cluster run, with the audit already used for the
   Context section above: `realtime`/`duration` ratio from `execution_trace.txt`. Expect the
   median slot to fall from ~79 s toward ~50 s (removing ~25 s of poll lag) and the useful-work
   fraction to rise from 19.3% toward ~35–40%. Throughput should rise from the observed 67
   tasks/min toward the low hundreds, bounded now by SLURM's ~96-per-QOS running cap and the
   ~5 tasks/s sbatch ceiling rather than by Nextflow's 100-slot window.
6. **`-resume` behaviour:** with `cache = 'lenient'` (`conf/crg.config:8`), removing
   `GUNZIP_TARGET_GFF` invalidates everything downstream of it — expect `FILTER_TARGET`,
   `GFF_STATS_TARGET` and their `GFFCOMPARE`s to re-run on the first resumed run. That is
   correct, not a regression.
