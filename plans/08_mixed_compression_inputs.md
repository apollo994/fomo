# Plan: Accept both compressed (.gz) and uncompressed FASTA / GFF3 inputs

## Context

The samplesheet schema (`assets/schema_input.json`) already *advertises* that
inputs may be gzipped or plain — both the `fasta` and `gff3` patterns end in
`(\.gz)?`. But the pipeline only works in practice when every input is gzipped:
four subworkflow steps call the nf-core `GUNZIP` module **unconditionally**, and
that module always runs `gzip -cd`, which aborts with "not in gzip format" on a
plain `.fasta`/`.gff3`. So a user who follows the schema and supplies an
uncompressed file gets a runtime crash. This plan closes the gap so a samplesheet
can freely mix compressed and uncompressed files per row.

## Current state (review findings)

**Unconditional GUNZIP — these 4 break on uncompressed input:**

| Site | File:line | What it decompresses |
|------|-----------|----------------------|
| `GUNZIP_FASTA` | `preprocessing.nf:24` | source FASTA |
| `GUNZIP_RAW_SOURCE_GFF` | `preprocessing.nf:126` | source GFF3 (for AGAT) |
| `GUNZIP_TARGET` | `projection.nf:19` | target FASTA (for minimap2) |
| `GUNZIP_TARGET_GFF` | `benchmarking.nf:15` | target GFF3 (for filter/AGAT) |

**Already compression-agnostic (no change needed):**
- `modules/local/filter_transcript.nf:20-24` — inline `if [[ $gff3 == *.gz ]]`
- `modules/local/gff_to_gene_bed.nf:22-26` — inline `if [[ $gff3 == *.gz ]]`
- `assets/schema_input.json` — already permits optional `.gz`

**Decompression cannot simply be dropped:** gffread, minimap2 (index/align),
samtools faidx all need an uncompressed FASTA; AGAT needs an uncompressed GFF3.
So a single conditional decompression point per input is the right design.

## Approach: a single `MAYBE_GUNZIP` local module (handles both inline)

Add `modules/local/maybe_gunzip.nf` — one process that takes `[meta, file]`,
decompresses it if it ends in `.gz`, otherwise copies it through, and in **both**
cases writes to the same `ext.prefix`-derived output name. Because the output
name is computed identically regardless of compression state, the intermediate
artefacts are byte-for-byte consistently named whether the user supplied a `.gz`
or a plain file — this is what eliminates the naming caveat.

This follows the same "handle `.gz` inline" pattern the repo already uses in
`modules/local/filter_transcript.nf` and `modules/local/gff_to_gene_bed.nf`,
rather than the nf-core GUNZIP module (which can only decompress and would need a
branch/mix wrapper plus a second rename step to normalise names).

```groovy
process MAYBE_GUNZIP {
    tag "${archive}"
    label 'process_single'

    // coreutils + gzip (same Wave image the nf-core gunzip module uses)
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(archive, stageAs: 'input/*')   // staged in a subdir so
                                                          // output never collides
                                                          // with an uncompressed input

    output:
    tuple val(meta), path("${output}"), emit: gunzip
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | head -1 | sed "s/^.*) //"'), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def isGz          = archive.name.endsWith('.gz')
    def nameWithoutGz = isGz ? archive.baseName : archive.name
    def extension     = file(nameWithoutGz).extension
    def name          = file(nameWithoutGz).baseName
    def prefix        = task.ext.prefix ?: name
    output = prefix + ".${extension}"
    """
    if [[ "${archive}" == *.gz ]]; then
        gzip -cd ${args} ${archive} > ${output}
    else
        cp -L ${archive} ${output}
    fi
    """

    stub:
    def isGz          = archive.name.endsWith('.gz')
    def nameWithoutGz = isGz ? archive.baseName : archive.name
    def extension     = file(nameWithoutGz).extension
    def name          = file(nameWithoutGz).baseName
    def prefix        = task.ext.prefix ?: name
    output = prefix + ".${extension}"
    """
    touch ${output}
    """
}
```

`stageAs: 'input/*'` puts the input in a subdirectory so the produced
`${output}` in the work-dir root can never collide with an uncompressed input of
the same name (the `cp` is always to a fresh path).

## Changes

### 1. New file
- `modules/local/maybe_gunzip.nf` — as above.

### 2. Swap the 4 includes (alias unchanged → call sites untouched)
In `preprocessing.nf`, `projection.nf`, `benchmarking.nf`, replace
`include { GUNZIP as <NAME> } from '../../modules/nf-core/gunzip/main'` with
`include { MAYBE_GUNZIP as <NAME> } from '../../modules/local/maybe_gunzip'`.

The four aliases — `GUNZIP_FASTA`, `GUNZIP_RAW_SOURCE_GFF`, `GUNZIP_TARGET`,
`GUNZIP_TARGET_GFF` — and their `<NAME>(...)` / `<NAME>.out.gunzip` usages stay
exactly as written (emit channel is still `gunzip`).

### 3. `conf/modules.config` — no selector path changes needed
Because `MAYBE_GUNZIP` is a **process** (not a subworkflow), the aliased process
path is identical to today's (`FOMO:PREPROCESSING:GUNZIP_RAW_SOURCE_GFF`, etc.).
The existing `withName` blocks and their `ext.prefix` closures keep applying
unchanged — and now apply to **both** compression states, which is exactly what
removes the naming caveat. No edits required here.

### 4. Fix stale doc comments
`projection.nf:12` and `benchmarking.nf:9` say `fasta.gz, gff3.gz`; change to
`fasta[.gz], gff3[.gz]`.

## Naming caveat — resolved by design

The old concern (pass-through keeping the original filename) does not arise:
`MAYBE_GUNZIP` always writes `${prefix}.${extension}` in both branches, so e.g.
`GUNZIP_RAW_SOURCE_GFF` yields `<id>.raw.gff3` whether the input was
`<id>.gff3.gz` or `<id>.gff3`. Intermediate names are now compression-invariant.

## Verification

1. **Compressed path unchanged:** `nextflow run main.nf -profile test,docker`
   (current `assets/samplesheet.csv`, all `.gz`) still completes; outputs identical.
2. **Uncompressed path now works:** create `assets/test_data` plain copies and a
   throwaway samplesheet, e.g.
   ```bash
   # decompress one source + the target in place (copies, keep .gz originals)
   for f in <target.fna.gz> <target.gff3.gz> <source.fna.gz> <source.gff3.gz>; do
       gunzip -k "$f"
   done
   ```
   Point a `samplesheet_plain.csv` at the decompressed files and run
   `nextflow run main.nf -profile test,docker --input assets/samplesheet_plain.csv`.
   Expect a clean run with the same module counts.
3. **Mixed path:** one samplesheet row compressed, another plain → both succeed in
   the same run (proves per-file branching, not a global toggle).
4. Confirm a `MAYBE_GUNZIP` task that received an uncompressed input produces an
   output named identically to the compressed case (e.g. `<id>.raw.gff3`), proving
   the prefix is applied in both branches.
5. `-resume` after switching a single file from `.gz` to plain re-runs only that
   input's `MAYBE_GUNZIP` task and its descendants.

## Out of scope

- bgzip vs plain gzip distinction (gunzip handles both; bgzip-specific indexing
  is not used on these inputs).
- Auto-recompressing outputs.
- Changing the schema — it already permits both forms.
