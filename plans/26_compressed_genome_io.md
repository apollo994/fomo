# Plan: keep genomes compressed — decompress on node-local disk, never on NFS

> **Status:** Implemented · **Date:** 2026-10-14
> **Decision:** Genomes stay compressed: gffread stages its genome on node-local scratch, there is no gunzip task and no `samtools faidx`.

**Status: implemented** — see **Verification** for what was checked here and what is
left for the full mollusca OctVul run.

## Problem

On the mollusca OctVul run (`cheesy_lagrange`, 2026-10-14, S=15, T=1,
`--include_mrna --include_decoy`), preprocessing took 84 min instead of 22 min on Oct 7,
and 40 tasks were killed at their time limit (exit 140), all in disk-heavy steps:
`EXTRACT_SEQUENCES` (20), `EXTRACT_DECOY_SEQUENCES` (13), `SAMTOOLS_FAIDX` (5) and
`GUNZIP_FASTA` (2). Some needed a third attempt.

The cause is how each uncompressed genome is read over NFS
(`/nfs/scratch01/.../work_fomo`):

- `GUNZIP_FASTA` writes the uncompressed genome into the shared work dir: 58 GB for 15
  genomes, from 17 GB of `.gz`.
- Every gffread task then reads that whole file **twice**. It first builds its own `.fai`,
  which means scanning the whole genome, then it reads the sequences. Measured on *Sepia
  lycidas* (5.2 GB): ~5 GB read per pass, ~45 MB/s on an idle NFS, under 1.5 s of CPU.
- Per species with both flags, that is 1 write and 9 full reads: `SAMTOOLS_FAIDX` ×1,
  `EXTRACT_SEQUENCES` ×2×2, `EXTRACT_DECOY_SEQUENCES` ×2×2. About 600 GB of NFS traffic
  per run, almost all of it in the first minutes of preprocessing. On a busy NFS day the
  passes crawl and hit the 15 min / 30 min limits.

A second, smaller source of the same pattern: `FILTER_TRANSCRIPT` decompresses each
source GFF3 into its work dir as `_input.gff` and never deletes it (~250 MB per task,
4.7 GB per OctVul run). `GFF_TO_GENE_BED` does the same.

## Measured alternative

Same Sepia genome and GFF3s, on a genoa64 compute node:

| | Time | NFS read |
|---|---|---|
| gffread on the NFS copy, `.fai` already present | 112–127 s | ~5 GB |
| as the pipeline runs it (builds the `.fai` first) | ~2× that; 14–30 min on Oct 14 | ~10 GB |
| `zcat` the `.gz` to node-local `/tmp`, gffread there (builds `.fai`) | 59 s + 10–30 s | 1.5 GB |

The extracted FASTA was byte-identical.

## Cluster facts this relies on (verified)

- Every SLURM job, batch jobs included, gets `TMPDIR=/tmp/$SLURM_JOB_ID`, set by
  `/etc/slurm/scripts/slurm.task_prolog`. `/tmp` is the node's local NVMe (872 GB) and
  is per-job (`JobContainerType`/`NamespaceType = namespace/tmpfs`), so SLURM removes it
  at job end, even after a SIGKILL.
- Nextflow forwards `TMPDIR` into the container (`SINGULARITYENV_TMPDIR` in
  `.command.run`), and CRG's Singularity has `mount tmp = yes`. Checked: inside the gffread
  image `$TMPDIR` exists and is writable, and `gzip`, `zcat`, `mktemp` and `awk` are present.
- Off the cluster (`docker`, plain local runs) `TMPDIR` may be unset. The code falls back
  to `/tmp`, the same convention Nextflow's own `.command.run` uses.

## Which tools read compressed input, and what changes

| Consumer | Input | Reads `.gz`? | Change |
|---|---|---|---|
| `MINIMAP2_INDEX` | target genome | yes (gzopen) | none |
| `GFF_STATS*` (gff-feature-stats) | GFF3 | yes | none |
| `CURATE_MODELS` | target reference GFF3 | yes (`gzip.open`) | none |
| `GFFREAD` (`EXTRACT_SEQUENCES`, `EXTRACT_DECOY_SEQUENCES`, `EXTRACT_PROJECTED_LNC`) | genome `-g` | **no** — random-access reader, no gzip/bgzf path | **decompress inside the task to `$TMPDIR`**, run gffread there |
| `SAMTOOLS_FAIDX` (chromosome sizes for the decoy chain) | genome | bgzip only; hard-errors on plain gzip (NCBI/Ensembl) | **replace** with `FASTA_SIZES`: one streaming pass, `gzip -cd \| awk` |
| `FILTER_TRANSCRIPT` | source GFF3 | no; the awk reads the file **twice** | decompress to `$TMPDIR`, not the work dir |
| `GFF_TO_GENE_BED` | source GFF3 | no; reads it once | stream: `gzip -cd \| awk … sizes -` under `pipefail` |

So the genomes stay `.gz` from the samplesheet until the moment a task needs plain text.
That plain copy lives on node-local disk for the task's lifetime and never touches NFS.

## Changes

1. **`modules/nf-core/gffread/main.nf` — FOMO-local deviation**, same pattern as
   `minimap2/align`. When `-g` is given, the genome is staged into a private
   `mktemp -d` directory under `${TMPDIR:-/tmp}`:
   - `.gz` → `gzip -cd`;
   - plain → `cp -L`, so NFS is read once instead of twice (gffread's `.fai` scan, then
     extraction).

   gffread runs against that copy and writes its `.fai` there; an `EXIT` trap removes
   the directory. `set -e` (Nextflow's default `bash -ue`) aborts the task with gzip's
   own exit code if decompression fails, so a corrupt archive cannot produce a truncated
   FASTA. Outputs, names and `ext.args`/`ext.prefix` are unchanged.

   The scratch directory falls back from `$TMPDIR`, to `/tmp`, to the task work dir.
   The last one is slower but still correct, and covers an executor whose `TMPDIR` is
   not visible inside the container. The `mktemp` template uses six `X`s: this image's
   `mktemp` rejected a four-`X` template when tested. Changes 3 and 4 share the same
   helper logic, inlined per module (no shared shell library exists in `bin/` for this).
2. **Delete `GUNZIP_FASTA` and `GUNZIP_TARGET`.** The gffread aliases take the samplesheet
   FASTA directly. `modules/local/maybe_gunzip.nf` then has no users and is deleted, along
   with its `FOMO:PROJECTION:GUNZIP_TARGET` selector.
3. **`FASTA_SIZES` replaces `SAMTOOLS_FAIDX`** in the decoy branch
   (`modules/local/fasta_sizes.nf`, container `ubuntu:22.04`, which has gzip and mawk).
   It emits `<name>\t<length>` in FASTA order. The name is the header up to the first
   whitespace, as `samtools faidx` does. Its output is byte-identical to
   `cut -f1,2 <genome>.fai`, the old `.sizes`. The now-unused
   `modules/nf-core/samtools/faidx` is deleted.
4. **`FILTER_TRANSCRIPT`**: `_input.gff` moves to a `mktemp -d` directory under
   `${TMPDIR:-/tmp}`, removed on exit. A plain input is still used in place (symlink).
5. **`GFF_TO_GENE_BED`**: no `_input.gff` at all. `gzip -cd` (or `cat`) is piped into
   the existing awk, which reads the sizes file and then `-`. The module already runs
   `set -euo pipefail`.
6. **Docs**: CLAUDE.md (task counts, the "no gunzip on the GFF3 side" paragraph, the decoy
   chain listing), the comments in `preprocessing.nf` / `projection.nf` that explained
   why the gunzips existed, and module comments that referred to "the MAYBE_GUNZIP trick".

### Effect (OctVul, both flags)

| | Before | After |
|---|---|---|
| uncompressed genome bytes written to NFS | 58 GB | 0 |
| genome bytes read from NFS | ~530 GB plain | ~85 GB compressed (5 reads × 17 GB) |
| `_input.gff` left in the work dir | ~4.7 GB | 0 |
| tasks | `S` `GUNZIP_FASTA` + `T` `GUNZIP_TARGET` + `S` `SAMTOOLS_FAIDX` | `S` `FASTA_SIZES` (decoy only) |
| CPU per gffread task | — | +25–60 s of `gzip -cd` (Sepia: 25 s user) |

## Not changed

- **Pipeline GFF3 outputs stay plain** (your choice). They are ~1.1 GB per OctVul run,
  and gzipping them would touch every consumer and the published file names.
- **`.mmi` indexes** (11 GB for T=1 here): minimap2 cannot read a compressed index, and
  dropping the separate index step was rejected in plans/24.
- **`MINIMAP2_ALIGN`'s intermediate SAM**: 0.2 GB total on OctVul, not worth touching.

## Risks

- **Local disk per task.** Each gffread task holds one uncompressed genome on the node's
  `/tmp`: 5 GB for the largest OctVul genome, ~15 GB for the 14.7 Gb *Vampyroteuthis*.
  With 826 GB free per node, about 50 concurrent 15 GB copies would fill one node. SLURM
  does not reserve `/tmp` space per job. If this ever bites, add `--tmp=` to the
  processes' `clusterOptions` so SLURM only places them on nodes with room.
- **Decompression is repeated** per gffread task: up to 4 per source with both flags,
  plus 1 for `FASTA_SIZES`. That is CPU traded for NFS bandwidth, which is the
  constraint here.
- **Resume**: the gffread, `FILTER_TRANSCRIPT` and `GFF_TO_GENE_BED` scripts change, so
  everything downstream of preprocessing and `EXTRACT_PROJECTED_LNC` re-runs on `-resume`.

## Verification

Done here:
- `FASTA_SIZES` on the Sepia `.gz` vs `cut -f1,2` of the gffread-built `.fai`: identical.
- `-profile test,crg --include_mrna --include_decoy` run with this change, against a
  baseline run of the previous commit with the same flags. Published outputs compared:
  `sources/*/annotation/*.filtered.gff3`, `targets/*/annotation/*.gff3`,
  `targets/*/curated/*`, `targets/*/gffcompare/*.stats`.
- `GUNZIP_*` and `SAMTOOLS_FAIDX` are absent from the new trace, and no task work dir
  holds an uncompressed genome or `_input.gff`.

For you on the OctVul set:
- preprocessing wall time, and no exit-140 rows on `EXTRACT_*` / `FASTA_SIZES`;
- `du -sh` of the work dir vs `cheesy_lagrange`.
