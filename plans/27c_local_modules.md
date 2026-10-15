# 27c — All modules local; containers pinned in one place; versions collected; no conda

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: U2, U4, U6.
Review items: C1–C5, D4 (FILTER_* merge), F15/F22 of the module review (versions, stubs).

## Status: not started

## Goal

The repository has one kind of module: local, trimmed to what FOMO calls, with one container
label, honest inputs, versions on the topic, and a working stub. `modules/nf-core/` and every
conda artefact are gone. Software versions show up in every MultiQC report.

## Steps

### A. Move and trim the nf-core modules

Each module becomes `modules/local/<process_lower>.nf`. Then update the includes and the config
selectors.

- [ ] **A1. `minimap2_index.nf` (MINIMAP2_INDEX).**
  - Input `tuple val(meta), path(fasta)`; output `tuple val(meta), path("*.mmi")`.
  - Add a tag. Keep the `-x splice -I 20G` args in config.
- [ ] **A2. `minimap2_align.nf` (MINIMAP2_ALIGN), trimmed from the fork.**
  - Input `tuple val(meta), path(query), path(mmi)`. Delete the four constant `val` inputs and
    `meta2`.
  - Script:
    - the exclude-self filter → `bin/fasta_drop_species.awk` (with a pytest);
    - `minimap2 ${args} -t … -a -o x.sam`;
    - `samtools sort --write-index` with CSI;
    - remove intermediates.
  - Outputs are BAM and CSI only. Remove the PAF, BAM-input and `samtools reset` branches.
  - `exclude_self` must apply on the only remaining code path.
  - Intermediates go to `${TMPDIR:-.}` scratch when available (plans/26 rule), *only if*
    `samtools sort` output stays identical (G4).
- [ ] **A3. `extract_transcripts.nf` (EXTRACT_TRANSCRIPTS), replacing the gffread fork.**
  - Input `tuple val(meta), path(gff), path(genome)`.
  - Script: stage the genome into node-local scratch (the same `mktemp -d` logic and trap), then
    `gffread -w ${prefix}.fasta -g genome gff`.
  - One output: `fasta`.
  - Aliases: `EXTRACT_SEQUENCES`, `EXTRACT_DECOY_SEQUENCES`, `EXTRACT_PROJECTED_LNC`. Keep them
    to limit churn, or rename to `EXTRACT_*_TRANSCRIPTS` if every selector is updated in the
    same step.
- [ ] **A4. `bedtools_complement.nf`, `samtools_stats.nf`, `seqkit_stats.nf`, `multiqc.nf`.**
  - Copy them, then drop unused inputs:
    - SAMTOOLS_STATS's `[[:],[],[]]` reference;
    - MULTIQC's logo, replace_names and sample_names.
  - Keep each tool's args exactly. The outputs MultiQC parses natively (samtools `.stats`) must
    not change.
- [ ] **A5. Delete `modules/nf-core/`** (including every `environment.yml`).

### B. Consolidate local modules

- [ ] **B1. Merge `FILTER_ANNOTATION` and `FILTER_GFF_BY_ID` into one module, `FILTER_GFF_BY_ID`.**
  - Input `tuple val(meta), path(gff, stageAs: 'input/*'), path(ids), path(fasta)`; `[]` when
    unused. The mode follows from which inputs are present (keep mode when an id map and FASTA
    are given).
  - Aliases: `FILTER_ANNOTATION` (keep) and `FILTER_ALLMODELS` (drop). The published
    `filtered.gff3` must be byte-identical (G4).
- [ ] **B2. Remove unused `when:` blocks.** Keep `when:` only on CURATE_MODELS, which is gated by
  `ext.when`. 27e revisits the gating.
- [ ] **B3. Remove unused inputs:** GFFCOMPARE_BATCH's `ref_meta`, and the `val(feature_type)` of
  FILTER_TRANSCRIPT, which duplicates `meta.feature_type`.
- [ ] **B4. Params are passed through `ext.args` only.** `RELOCATE_LOCI` stops reading
  `params.relocate_seed` and `params.decoy_cap` inside the module; config builds the args.

### C. Containers (U6, T5)

- [ ] **C1. Write `conf/containers.config`, included after `base.config`.** One block
  `withLabel: 'img_<name>' { container = { workflow.containerEngine in ['singularity','apptainer'] ? '<singularity ref>' : '<docker ref>' } }`
  per image family:

  | label | image | modules |
  |:---|:---|:---|
  | `img_python` | `python:3.11@sha256:…` | python stdlib scripts |
  | `img_ubuntu` | `ubuntu:22.04@sha256:…` | awk / coreutils |
  | `img_samtools` | Wave `htslib_samtools:1.23.1` | |
  | `img_minimap2` | Wave `minimap2:2.30` (index) | |
  | `img_minimap2_samtools` | | align |
  | `img_gffread` | | |
  | `img_bedtools` | | |
  | `img_seqkit` | | |
  | `img_multiqc` | biocontainers 1.35 | |
  | `img_curate` | Wave `python_gffcompare_gffread` | |
  | `img_td2` | Wave td2 1.1.0 | |
  | `img_gffcompare` | biocontainers 0.12.6 | until 27g moves GFFCOMPARE_BATCH to `img_curate` |

  Keep the exact images in use today; only pin the floating tags.
- [ ] **C2. Pin the floating tags by digest.** Get them from the registry API, e.g.:

  ```
  curl -s "https://auth.docker.io/token?service=registry.docker.io&scope=repository:library/python:pull" | jq -r .token
  ```

  then `HEAD /v2/library/python/manifests/3.11` with
  `Accept: application/vnd.oci.image.index.v1+json`. Record the digest and the date in a
  comment in `containers.config`. For Singularity use `docker://python:3.11@sha256:…`.
- [ ] **C3. Each module carries exactly one `img_*` label, and no `container` or `conda`
  line.** Add `tests/test_repo_conventions.py`, which greps the modules and fails otherwise.
- [ ] **C4. The `docker` profile's amd64 pin selects `withLabel: 'amd64_only'`.** Add that label
  to every module that runs the vendored `gff-feature-stats`. Delete the name regex.
- [ ] **C5. Drop the `conda` profile and `conda.channels`** from `nextflow.config`, and every
  `conda` directive.
- [ ] **C6. Pull every pinned image once with Singularity** (cache `/software/rg/images`). G4 with
  the pinned digests must be identical; the digest of `python:3.11` may differ from the image
  used for the baseline. If G4 changes, stop: a library difference inside the image matters.

### D. Versions (U4)

- [ ] **D1. Every module emits `topic: versions`, and no eval fails.**
  - Add versions to filter_gff_by_id (merged), gff_by_source, merge_fasta, seqkit_to_mqc
    (deleted in 27i; skip if not worth it), and td2_coding_filter.
  - Fix the awk evals: `awk -W version 2>&1 | head -1` for mawk, or report the base image
    version. Verify inside `img_ubuntu`.
  - FILTER_TRANSCRIPT moves from `versions.yml` to the topic.
- [ ] **D2. Collect the versions in `fomo.nf`.**
  - `channel.topic('versions')` → unique (process, tool, version) → YAML in MultiQC's
    `software_versions` format, through `collectFile` as `software_versions_mqc.yml`.
  - Pass the file to every per-target MULTIQC and to MULTIQC_RUN_SUMMARY.
  - Put a helper function in `subworkflows/local/fomo_utils.nf` (created here if 27e has not
    created it yet).
- [ ] **D3. MULTIQC's own version goes on the topic as well.**
- [ ] **D4. Confirm the manifest floor `!>=24.10.0` supports `channel.topic` without a preview
  flag.** Check the Nextflow changelog or docs. If it does not, raise the floor in 27d and
  record it here.

### E. Tags, stubs and the vendored binary

- [ ] **E1. Every module has a stub, and S1–S3 still pass.** Use one tag rule (master §8); 27e
  refines it once `meta.track` exists.
- [ ] **E2. Vendored binary provenance.**
  - Add `bin/gff-feature-stats.sha256` and `bin/gff-feature-stats.PROVENANCE` (upstream repo,
    commit `67ffea7`, build command, glibc floor).
  - pytest checks the checksum.
  - Remove the provenance prose from module comments.

## Acceptance

- G1, G2, G3 and G4 identical.
- `modules/nf-core` does not exist; `test_repo_conventions.py` passes.
- `grep -rn "conda" modules conf nextflow.config` finds no directive.
- V1's MultiQC data contains a software-versions table that lists every tool. This is a new
  section, which U1 allows.
- Results list the digests and dates, and the module line count before and after.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|

## Results

*(fill in)*
