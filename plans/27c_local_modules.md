# 27c — All modules local; containers pinned in one place; versions collected; no conda

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: U2, U4, U6.
Review items: C1–C5, D4 (FILTER_* merge), F15/F22 of the module review (versions, stubs).

## Status: done (2026-10-15)

## Goal

The repository has one kind of module: local, trimmed to what FOMO calls, with one container
label, honest inputs, versions on the topic, and a working stub. `modules/nf-core/` and every
conda artefact are gone. Software versions show up in every MultiQC report.

## Steps

### A. Move and trim the nf-core modules

Each module becomes `modules/local/<process_lower>.nf`. Then update the includes and the config
selectors.

- [x] **A1. `minimap2_index.nf` (MINIMAP2_INDEX).**
  - Input `tuple val(meta), path(fasta)`; output `tuple val(meta), path("*.mmi")`.
  - Add a tag. Keep the `-x splice -I 20G` args in config.
- [x] **A2. `minimap2_align.nf` (MINIMAP2_ALIGN), trimmed from the fork.** (updated by D27c-1,
  D27c-2)
  - Input `tuple val(meta), path(query), path(mmi)`. Delete the four constant `val` inputs and
    `meta2`.
  - Script:
    - the exclude-self filter → `bin/fasta_drop_species.sh` (with a pytest; D27c-1);
    - `minimap2 ${args} -t … -a -o x.sam`;
    - `samtools sort --write-index` with CSI;
    - remove intermediates.
  - Outputs are BAM and CSI only. Remove the PAF, BAM-input and `samtools reset` branches.
  - `exclude_self` must apply on the only remaining code path.
  - Intermediates go to `${TMPDIR:-.}` scratch when available (plans/26 rule), *only if*
    `samtools sort` output stays identical (G4).
- [x] **A3. `extract_transcripts.nf` (EXTRACT_TRANSCRIPTS), replacing the gffread fork.**
  - Input `tuple val(meta), path(gff), path(genome)`.
  - Script: stage the genome into node-local scratch (the same `mktemp -d` logic and trap), then
    `gffread -w ${prefix}.fasta -g genome gff`.
  - One output: `fasta`.
  - Aliases: `EXTRACT_SEQUENCES`, `EXTRACT_DECOY_SEQUENCES`, `EXTRACT_PROJECTED_LNC`. Keep them
    to limit churn, or rename to `EXTRACT_*_TRANSCRIPTS` if every selector is updated in the
    same step.
- [x] **A4. `bedtools_complement.nf`, `samtools_stats.nf`, `seqkit_stats.nf`, `multiqc.nf`.**
  - Copy them, then drop unused inputs:
    - SAMTOOLS_STATS's `[[:],[],[]]` reference;
    - MULTIQC's logo, replace_names and sample_names.
  - Keep each tool's args exactly. The outputs MultiQC parses natively (samtools `.stats`) must
    not change.
- [x] **A5. Delete `modules/nf-core/`** (including every `environment.yml`). The forked gffread
  and minimap2/align comments were not swept in 27b (D27b-3): write fresh, policy-conform
  comments for the local versions.

### B. Consolidate local modules

- [x] **B1. Merge `FILTER_ANNOTATION` and `FILTER_GFF_BY_ID` into one module, `FILTER_GFF_BY_ID`.**
  - Input `tuple val(meta), path(gff, stageAs: 'input/*'), path(ids), path(fasta)`; `[]` when
    unused. The mode follows from which inputs are present (keep mode when an id map and FASTA
    are given).
  - Aliases: `FILTER_ANNOTATION` (keep) and `FILTER_ALLMODELS` (drop). The published
    `filtered.gff3` must be byte-identical (G4).
- [x] **B2. Remove unused `when:` blocks.** Keep `when:` only on CURATE_MODELS, which is gated by
  `ext.when`. 27e revisits the gating.
- [x] **B3. Remove unused inputs:** GFFCOMPARE_BATCH's `ref_meta`, and the `val(feature_type)` of
  FILTER_TRANSCRIPT, which duplicates `meta.feature_type`.
- [x] **B4. Params are passed through `ext.args` only.** `RELOCATE_LOCI` stops reading
  `params.relocate_seed` and `params.decoy_cap` inside the module; config builds the args.

### C. Containers (U6, T5)

- [x] **C1. Write `conf/containers.config`, included after `base.config`.** (updated by D27c-3) One block
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
- [x] **C2. Pin the floating tags by digest.** (updated by D27c-3) Get them from the registry API, e.g.:

  ```
  curl -s "https://auth.docker.io/token?service=registry.docker.io&scope=repository:library/python:pull" | jq -r .token
  ```

  then `HEAD /v2/library/python/manifests/3.11` with
  `Accept: application/vnd.oci.image.index.v1+json`. Record the digest and the date in a
  comment in `containers.config`. One reference for both engines, `python@sha256:…` (no tag:
  Singularity rejects tag + digest, D27c-3).
- [x] **C3. Each module carries exactly one `img_*` label, and no `container` or `conda`
  line.** Add `tests/test_repo_conventions.py`, which greps the modules and fails otherwise.
- [x] **C4. The `docker` profile's amd64 pin selects `withLabel: 'amd64_only'`.** Add that label
  to every module that runs the vendored `gff-feature-stats`. Delete the name regex.
- [x] **C5. Drop the `conda` profile and `conda.channels`** from `nextflow.config`, and every
  `conda` directive.
- [x] **C6. Pull every pinned image once with Singularity** (cache `/software/rg/images`). G4 with
  the pinned digests must be identical; the digest of `python:3.11` may differ from the image
  used for the baseline. If G4 changes, stop: a library difference inside the image matters.

### D. Versions (U4)

- [x] **D1. Every module emits `topic: versions`, and no eval fails.**
  - Add versions to filter_gff_by_id (merged), gff_by_source, merge_fasta, seqkit_to_mqc
    (deleted in 27i; skip if not worth it), and td2_coding_filter.
  - Fix the awk evals: `awk -W version 2>&1 | head -1` for mawk, or report the base image
    version. Verify inside `img_ubuntu`.
  - FILTER_TRANSCRIPT moves from `versions.yml` to the topic.
- [x] **D2. Collect the versions in `fomo.nf`.** (updated by D27c-4, D27c-5)
  - `channel.topic('versions')` → unique (process, tool, version) → YAML in MultiQC's
    `software_versions` format, through `collectFile` as `fomo_software_mqc_versions.yml`
    (D27c-4).
  - Pass the file to every per-target MULTIQC and to MULTIQC_RUN_SUMMARY, and to nothing else
    (D27c-5).
  - Put a helper function in `subworkflows/local/fomo_utils.nf` (created here if 27e has not
    created it yet).
- [x] ~~**D3. MULTIQC's own version goes on the topic as well.**~~ Not possible (D27c-5): MULTIQC
  emits no versions; the report shows its own MultiQC version.
- [x] **D4. Confirm the manifest floor `!>=24.10.0` supports `channel.topic` without a preview
  flag.** Check the Nextflow changelog or docs. If it does not, raise the floor in 27d and
  record it here. → It does not (D27c-6): 27d raises the floor to `!>=25.04.0`.

### E. Tags, stubs and the vendored binary

- [x] **E1. Every module has a stub, and S1–S3 still pass.** Use one tag rule (master §8); 27e
  refines it once `meta.track` exists. (updated by D27c-7: the tag rule moves to 27e B2.)
- [x] **E2. Vendored binary provenance.**
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
| D27c-1 | A2 | exclude-self filter in `bin/fasta_drop_species.awk` | `bin/fasta_drop_species.sh` (bash, `set -euo pipefail`, usage check, then the same awk one-liner) + `tests/fasta_drop_species/` | an awk shebang needs a fixed awk path in every image; a bash script follows the master §8 shell convention and is shellcheck-linted | none |
| D27c-2 | A2 | intermediates to `${TMPDIR}` scratch if output stays identical | not done: the unsorted SAM and the filtered query stay in the work dir, as before | minimap2 and samtools write their command lines into the BAM `@PG` header, so a scratch path would put a random `mktemp` name into every published BAM header; avoiding it needs stdin/stdout redirection and a comment explaining why. The SAM is written once and read once (a few hundred MB), unlike the genome that plans/26 moved off NFS. The step was conditional. | none |
| D27c-3 | C1, C2 | Singularity uses `docker://python:3.11@sha256:…` | `python@sha256:…` and `ubuntu@sha256:…` (tag in a comment), one string for both engines | `singularity pull docker://ubuntu:22.04@sha256:…` → "Docker references with both a tag and digest are currently not supported" (tested 2026-10-15); digest-only pulls fine and is glibc 2.35 | none |
| D27c-4 | D2 | YAML named `software_versions_mqc.yml` | `fomo_software_mqc_versions.yml`, `{process: {tool: version}}` plus a `Workflow` entry (Nextflow, fomo) | MultiQC 1.35 finds versions files by `fn_re: .+_mqc_versions\.(yaml|yml)`; a `*_mqc.yml` name is read as custom content. Tested in the 1.35 image: renders the Software Versions section and `multiqc_software_versions.txt` | 27i B3 passes this file through unchanged |
| D27c-5 | D2, D3 | MULTIQC's version on the topic; YAML passed to the reports | MULTIQC emits no versions. The YAML reaches the MULTIQC tasks only (mixed into the per-target broadcast files in REPORTING and into MULTIQC_RUN_SUMMARY's inputs), never `ch_all_mqc` | a process that waits on the versions topic and also emits to it never lets the topic close (deadlock). `ch_all_mqc` feeds RUN_SUMMARY_TABLES, which emits python's version | 27g/27h/27i: the versions YAML must never be an input of a process that emits to the topic (RUN_SUMMARY_TABLES, any renderer outside MULTIQC) |
| D27c-6 | D4 | confirm `!>=24.10.0` supports topics | it does not | topic channels needed `nextflow.preview.topic` in 24.04 and 24.10 and left preview in 25.04 (Nextflow 25.04 migration notes, docs.seqera.io/nextflow/migrations/25-04; topic-channels tutorial). The pre-27c modules already used `topic: versions`, so 24.10 was never a working floor | 27d A7 raises `nextflowVersion` to `!>=25.04.0` (updated by D27c-6) |
| D27c-7 | E1 | one tag rule in every module | new modules use the existing `<target>.<feature_type>[.decoy]` / `<id>.<feature_type>[.decoy]` forms; the unified rule is applied once, in 27e B2 | the master §8 rule needs `meta.track`, which 27e introduces; applying an interim rule now and the final one in 27e would touch every module twice | 27e B2 applies the rule to every module (no change to its scope) |

## Results

Gates (2026-10-15):

- G1 `pytest -q`: 122 passed (87 + 4 `fasta_drop_species` + 31 `test_repo_conventions`).
- G2 `pre-commit run --all-files`: ruff, shellcheck, `nextflow lint` (43 files, 0 errors), trailing
  whitespace — all pass.
- G3 `nf-test test tests/pipeline --profile +singularity`: S1–S4 pass. The S1 work dir holds a
  `fomo_software_mqc_versions.yml` with every tool and no failed eval (mawk `1.3.4 20200120`,
  coreutils `8.32`, gff-feature-stats `0.2.0`, …).
- G4 `tests/golden/gate.sh …/gate/27c`: V1, V2, V3, V4 identical — with the newly pinned
  python (3.11.17) and ubuntu digests, so C6 holds. Per-process cpus/memory/time in the V1 trace
  are identical to the 27b gate run.
- MultiQC: both V1 target reports and the run summary have a Software Versions section
  (`multiqc_software_versions.txt`) listing every process and the `Workflow` row.
- Cluster use: 0.9 CPU-h (nf-test 9 min, mostly image pulls; G4 4 × 1 min; 4 CPUs each).

Modules: before, 28 `.nf` files (21 local + 7 nf-core) plus 7 `environment.yml`, 1,272 lines;
after, 27 `.nf` files in `modules/local`, 964 lines. `modules/nf-core` is gone; `grep -rn conda
modules conf nextflow.config` finds nothing.

Images (`conf/containers.config`, one block per `img_*` label; 12 labels):

| label | reference | note |
|:---|:---|:---|
| `img_python` | `python@sha256:50cfccc707e8721201ee018bfd7f6a85d8eacf49c4f9f75833a2aba09a10411f` | `python:3.11` index digest, 2026-10-15 (Python 3.11.17) |
| `img_ubuntu` | `ubuntu@sha256:5ec03bb3441e8b0bf3b4f9cd4629a1ae763010dc3035bb8da3ae6cf026486401` | `ubuntu:22.04` index digest, 2026-10-15 (22.04.5, glibc 2.35) |
| others | the Wave / BioContainers references in use before 27c | unchanged |

Moved images: GFFCOMPARE_TO_SCATTER and SAMTOOLS_TO_MQC from the MultiQC image to `img_python`
(stdlib scripts); MERGE_FASTA from the Wave coreutils image to `img_ubuntu` (that image is no
longer used). Versions were added to FILTER_GFF_BY_ID, GFF_BY_SOURCE, MERGE_FASTA, SEQKIT_TO_MQC,
TD2_CODING_FILTER and FILTER_TRANSCRIPT (topic instead of `versions.yml`); the three broken
`awk --version` evals now report mawk.

Other: `bin/gff-feature-stats.sha256` + `bin/gff-feature-stats.PROVENANCE` (checksum tested);
FILTER_ANNOTATION is an alias of FILTER_GFF_BY_ID (keep mode); GFFCOMPARE_BATCH, MINIMAP2_ALIGN,
EXTRACT_TRANSCRIPTS, BEDTOOLS_COMPLEMENT and SAMTOOLS_STATS take one tuple each, so the four
`multiMap` call sites are plain `map`s. V1 logs one more "no process matching config selector"
warning (RELOCATE_LOCI's new `ext.args` block, decoy-only): the D27a-4 class, handled in 27e C2.
