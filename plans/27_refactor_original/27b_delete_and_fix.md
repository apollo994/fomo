# 27b — Delete dead code, fix latent bugs, comment sweep

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Review items: §2 (code
ablation), R1–R3, P1 (species pattern), V1–V4, D6.

## Status: not started

## Goal

Remove everything that is dead, unused or redundant, which should change no deliverable. Fix the
latent bugs the review confirmed. Bring every comment and user-facing text in line with the
comment policy (master §8), so later subplans start from present-tense, accurate code.

## Steps

### A. Delete

- [ ] **A1. `modules/nf-core/gunzip` and `modules/nf-core/gffcompare`.**
  - Remove the hidden dependencies first:
    - the `gffcompare_batch.nf` conda directive pointing into `../nf-core/gffcompare` (conda is
      dropped in 27c anyway; delete the line now);
    - `merge_fasta.nf` borrows the gunzip image. Keep the image string; containers move in 27c.
  - Fix the `curate_models.nf` comment that names these modules.
- [ ] **A2. Tracked junk and unused assets:**
  - `modules/nf-core/.DS_Store`;
  - `assets/samplesheet_full.csv`;
  - `assets/test_data/Lycaena_phlaeas_282391/` and `Lycaena_virgaureae_282395/` (in no
    samplesheet);
  - update `bin/subsample_test_data.sh`'s `PAIRS` list to the three species kept.
- [ ] **A3. Dead config:**
  - `withLabel: 'big_mem'/'small'` and the commented GPU block in `crg.config`;
  - the MINIMAP2_INDEX `ext.prefix` (never read);
  - the selector alternatives `GFFCOMPARE|` (`modules.config` publish rule) and `_PROJECTED`
    without `_BATCH` (docker profile);
  - the six `ext.prefix` blocks that restate the module default (FILTER_TRANSCRIPT,
    GFF_TO_GENE_BED, RELOCATE_LOCI, FILTER_ANNOTATION, SEQKIT_TO_MQC, CURATE_MODELS). Verify
    each default really is identical first;
  - the `versions.yml` `saveAs` filters in the 8 publishDir blocks. No published process writes
    one; 27c moves FILTER_TRANSCRIPT to the topic.
- [ ] **A4. Unused emits and their cardinality comments:**
  - `PREPROCESSING.intergenic_bed`, `decoy_gff3`;
  - `PROJECTION.bam`, `index`;
  - `CURATION.curated`, `merged`, `report`;
  - the REPORTING and RUN_SUMMARY emits;
  - the unused `TD2_NONCODING.coding_ids` consumer in preprocessing. The subworkflow emit stays,
    because projection uses it.

  Check that every one is unused, with grep and a G3 run.
- [ ] **A5. Per-source re-key in `benchmarking.nf:77-86`:**
  - replace the block with `GFFCOMPARE_BATCH.out.stats.transpose()` (the reconstructed `id` is
    read by nothing);
  - delete the three comments about the "template and arithmetic must agree" coupling
    (`modules.config` around the split template, `benchmarking.nf`, `gff_stats_projected_batch.nf`).
- [ ] **A6. String-key join in `projection.nf:198-206`:** join `ch_allmodels_raw` with
  `SPLIT_GFF_BY_SOURCE.out.gff3` on meta, with `remainder: true`. Verify that both sides carry
  the *same* meta map (equal keys and values). The S1–S3 stub tests and V3 cover the empty path.
- [ ] **A7. Dead Python:**
  - `relocate_loci.build_shape`;
  - `Feature.attrs_raw`;
  - `fomo_stats.Identity.is_aggregate` and `is_self`. Either delete them, or replace the inline
    re-implementations in `run_summary_tables.py` with them; choose the smaller diff;
  - the discarded second return of `curate_models.read_models`;
  - the `top\d+` alternative in the filename inference;
  - `gff_stats_to_mqc --section-id`, which is only logged. Drop the flag from the module call
    too.
- [ ] **A8. `assets/empty_ids.txt`:** make the drop list of `FILTER_GFF_BY_ID` optional
  (`path(drop_ids)`, with `[]` for non-lncRNA tracks), and have `filter_gff_by_id.py` treat a
  missing `--drop` as "drop nothing". Delete the asset.

### B. Fix latent bugs (each with a pytest that failed before the fix)

- [ ] **B1. (R1) Batch modules mask failures.**
  - In `gffcompare_batch.nf` and `gff_stats_projected_batch.nf`, the worker function starts with
    `set -euo pipefail`.
  - Add exit codes 123 and 125 to the retry list *for these two processes only* in
    `base.config`.
  - Test: a shell test that runs the generated loop pattern with a failing first command, and
    expects a non-zero exit.
- [ ] **B2. (R2) `run_summary_tables.curation_index`.**
  - Merge the observed and novel rows into separate namespaces, or into two dicts.
  - Make the docstring true.
  - Test: a fixture where observed `genes` ≠ novel `genes` must report the novel value.
    `run_summary.json` stays identical (G4).
- [ ] **B3. (R3) `bam_to_gff.sh`.** A read name that does not have exactly three `|` fields
  must exit 1 with the read name in the message. Test with a fixture BAM.
- [ ] **B4. (U5) Species pattern.**
  - In `assets/schema_input.json`, set `species` to `^[A-Za-z0-9_-]+$` with a clear
    `errorMessage`.
  - Add `"uniqueEntries": ["species"]` to the array schema, and delete the hand-written
    duplicate check in `fomo.nf` if nf-schema 2.3 reports duplicates clearly. Otherwise keep it
    and record why.
  - Test: an nf-test pipeline test that expects a launch failure for a sheet with a species
    `Sp.1`.
- [ ] **B5. If 27a recorded a V4 ≠ V1 difference, fix it here** (it is a real bug: plain gzip
  vs bgzip).

### C. Comment and text sweep (master §8)

- [ ] **C1. User-facing text.**
  - `nextflow_schema.json` help:
    - `input` drops "top-N consensus";
    - `curate` describes one transcript per chain (t1 = most species → splice junctions →
      length) and has no "See plans/21";
    - `td2_psauron_min` says "pooled per target".
  - `assets/multiqc/sections.yml`:
    - the samtools_align description says per (target, track), not per source;
    - the td2_coding description has no "per-source" and no "consensus";
    - the header lists all sections.
  - `main.yml` and `run_summary_main.yml`: "at least two species" becomes "at least
    `--curate_min_species` species".
- [ ] **C2. Stale comments that contradict the code.** Fix or delete every one listed in review
  §3.5 V1, including:
  - `fomo.nf:83-84`;
  - `preprocessing.nf:20`;
  - `modules.config:110-111`;
  - the `benchmarking.nf` "2 aggregate models";
  - the `reporting.nf` consensus mentions;
  - `td2_noncoding.nf:7`;
  - `gff_stats_projected_batch.nf:50-51`;
  - `gffcompare_batch.nf:27`;
  - `gffcompare_accuracy_mqc.py` (FILTER_TARGET);
  - `td2_predict.nf:70-74` ("norm for self-projection");
  - `nextflow.config:51-53`;
  - `nextflow.config:126-127` (until 27d implements timestamps, say "overwritten each run");
  - `crg.config:11,18-19`;
  - `test.config:7,16`;
  - every `file:line` reference.
- [ ] **C3. Run and incident references (V2).** Delete every reference to runs, job IDs, clades,
  measurements and incidents. Where a rule is worth keeping, keep the rule in one line. Check
  with:

  ```
  grep -rniE 'jobid|job [0-9]{6,}|agitated_|thirsty_|mollusca|lepidoptera|apiaceae|cephalopod|vampyroteuthis|tick chromosome|seen for real|verified the hard way|82 targets' -- workflows subworkflows modules conf nextflow.config bin
  ```

  The command must return nothing.
- [ ] **C4. History narration (V3).** Rewrite in present tense or delete. Check with:

  ```
  grep -rniE 'used to|no longer|the old |un-batched|since plans|pre-plans|went with|is gone|are gone' -- <same paths>
  ```

  The command must return nothing.
- [ ] **C5. Verbosity (V4).**
  - Cut comments longer than about 3 lines that do not state an invariant.
  - Remove comments that describe absent things (`modules.config:279-283,286-288`).
  - Reduce duplicate explanations (D5) to one canonical comment. CLAUDE.md is rewritten in 27j;
    until then it may still repeat them.
  - Record the comment ratio before and after per layer
    (`python3 ablation/scripts/code_metrics.py`).

## Acceptance

- G1, G2, G3 (the stub matrix may lose expected emits but not tasks), and G4 identical.
- The C3 and C4 greps return nothing.
- The new tests for B1–B4 fail on the pre-fix code (show it once) and pass after.
- Results record the lines deleted (`git diff --stat` against the 27a commit) and the comment
  ratio per layer.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|

## Results

*(fill in)*
