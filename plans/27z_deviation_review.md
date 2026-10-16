# 27z — Deviation review of the 27 refactor

Part of [27_refactor.md](27_refactor.md); written in [27k](27k_final_validation.md) (D1–D4).

## Status: done (2026-10-16)

## Method

- **D1.** `diff -u` of every `plans/27_refactor_original/27*.md` against the working copy, saved to
  `/nfs/scratch01/rg/fzanarello/fomo_refactor/plan_diff.txt` (1,697 lines, 65 hunks). (`diff -ru
  --include` does not exist in GNU diff, so the files are diffed one by one.)
- **D2.** Every hunk is classified as bookkeeping (status, checkboxes, Results, status log), a
  documented deviation (the changed step cites a `D27x-N` row) or an undocumented change.
- **D3.** Every deviation is checked against U1–U8, the invariants (master §4) and the gate
  (master §5), and its evidence is located. A deviation is valid only if all three hold.
- **D4.** Every step of every original subplan 27a–27j is checked against the repository at
  `67ccc93` (not against the checkboxes): done, replaced or cut by a deviation, or superseded by a
  later subplan.

## D2. Hunk classification

| file | hunks | bookkeeping | documented deviation | undocumented → fixed in 27k |
|:---|---:|---:|---:|:---|
| 27_refactor.md | 3 | 3 | 0 | T2 still said GFFCOMPARE_BATCH moves to `img_curate` and that the record emitters are awk scripts; invariant 8 still named 27a as the owner of the TD2 test → notes citing D27g-1, D27g-3 and D27j-3 |
| 27a | 6 | 6 | 0 (4 rows, findings and tooling) | — |
| 27b | 10 | 10 | 0 (4 rows) | — (the A2 slip is D27f-1, with a correction note in Results) |
| 27c | 7 | — | A2, A5, C1, C2, D2, D3, D4, E1 | — |
| 27d | 7 | — | A5, A6, A7 | — |
| 27e | 4 | — | A1, B2, C1, C2 | C3 was struck through as "moot" without an id → **D27e-3** (retroactive) |
| 27f | 5 | — | B3 | — |
| 27g | 4 | — | A2, B2, C2, D1, D2 | the design section above the steps (identity built in each module script; `.awk` emitters; a Python gffcompare driver) did not show D27g-1/3/5 → notes added |
| 27h | 4 | — | — | A1 carried "(updated by D27g-5)" for a fact from 27g's Results (where the fixture records live), not a deviation → relabelled |
| 27i | 6 | — | B2, B3, C2, C4 | C7 still said "keep `asList`", which D27e-1 never created → note citing D27e-1 |
| 27j | 7 | — | C2, D3, D4 | — |
| 27k | 2 | 1 | B2 (D27g-3) | — |

Every hunk that changes the text of a step cites a deviation; the rest are checkboxes, statuses,
Results and the status log. The four undocumented changes above are propagation misses of
deviations that were themselves recorded; the only change of substance without an id was 27e C3,
now D27e-3. Nothing was reverted: none of them was a mistake in the code.

## D3. Deviations

Classes: **correction** (the plan was wrong or impossible as written), **improvement** (a better
route to the same goal), **scope cut** (a step not done), **scope addition** (work the plan did not
ask for).

| id | subplan | summary | class | consistent with U1–U8? | invariants and gate preserved? | evidence present? | verdict |
|:---|:---|:---|:---|:---|:---|:---|:---|
| D27a-1 | 27a | `tests/golden/gate.sh` submits V1–V4 and compares all four | scope addition | yes (U3) | yes: it runs exactly G4 | `tests/golden/gate.sh`; every subplan's Results | valid |
| D27a-2 | 27a | `ruff format` configured early, on new test files only; py311 lint codes in the debt list | improvement | yes (U3) | yes | `pyproject.toml`, `.pre-commit-config.yaml`; 27f removed the debt | valid |
| D27a-3 | 27a | stub path lists are a subset of the real path set (stubs cannot know source names) | correction | yes | yes: the real path set is G4's | `tests/pipeline/expected/*.paths` | valid |
| D27a-4 | 27a | finding: `if`-gated included processes warn "no process matching config selector" | correction (finding) | yes | — | 27a/27d V1 logs, quoted in D27e-2 | valid |
| D27b-1 | 27b | `Identity.is_aggregate`/`is_self` deleted (one of the two options the step offered) | scope cut | yes | yes | `git show 6b65e4b`; `Identity` itself went in 27i | valid |
| D27b-2 | 27b | the B2 test cannot fail before the fix: the bug was latent | correction | yes | yes | `test_curated_genes_come_from_the_novel_table`; B1/B3/B4 failed before (row lists the outputs) | valid |
| D27b-3 | 27b | forked nf-core module comments left to 27c, which rewrote the modules | scope cut (moved) | yes | yes | 27c A5 note; the modules are gone | valid |
| D27b-4 | 27b | the review's regexes flag present-tense phrases; the budget uses the narrower 27b lists | correction | yes | yes | `dev/code_metrics.py` `BUDGET_*`; the review regexes find 5 present-tense hits at `67ccc93` (27k C1) | valid |
| D27c-1 | 27c | `fasta_drop_species.sh` (bash + awk) instead of an awk script | correction | yes | yes (invariant 3: G4, stub matrix) | `bin/fasta_drop_species.sh`, `tests/fasta_drop_species/` | valid |
| D27c-2 | 27c | MINIMAP2_ALIGN intermediates stay in the work dir (conditional step) | scope cut | yes | yes: a scratch path would change the BAM `@PG` header | reasoning in the row; step was conditional | valid |
| D27c-3 | 27c | digest-only image references (`python@sha256:…`) for both engines | correction | yes (U6: pinned) | yes (G4 identical with the pinned digests) | the Singularity error quoted in the row; `conf/containers.config` | valid |
| D27c-4 | 27c | versions file `fomo_software_mqc_versions.yml` (MultiQC's `_mqc_versions` pattern) | correction | yes (U4) | yes | Software Versions section in every report (27c, 27i Results) | valid |
| D27c-5 | 27c | MULTIQC emits no versions; the YAML reaches only MULTIQC tasks (topic deadlock) | correction | yes (U4: the section exists; MultiQC shows its own version) | yes | the deadlock reasoning; G3/G4 complete | valid |
| D27c-6 | 27c | the `!>=24.10.0` floor never worked with topics → 25.04 | correction | yes: the code already needed 25.04, the manifest was wrong; not a behaviour change under U5 | yes | Nextflow 25.04 migration notes; CHANGELOG states it | valid |
| D27c-7 | 27c | the unified tag rule moved to 27e B2 | scope cut (moved) | yes | yes | 27e B2 done; `test_module_shape` | valid |
| D27d-1 | 27d | timestamp via hidden param `trace_report_suffix` (strict config syntax) | correction | yes (it implements U5) | yes | `nextflow.config`, schema, `docs/output.md`; two runs keep two report sets (27d Results) | valid |
| D27d-2 | 27d | `manifest.contributors` instead of the deprecated `author` | correction | yes | yes | `nextflow.config` | valid |
| D27d-3 | 27d | TD2 told apart by `meta.stage`; the coding-id join keyed on `target_id` | correction | yes | yes (G4 identical) | `projection.nf`; 27d Results | valid |
| D27d-4 | 27d | per-process comparison on effective values from a stub `crg` run | improvement | yes | yes: stronger than the flat dump | `fomo_refactor/config_{before,after}/stub_crg/resolved.tsv`, `resolve.sh` | valid |
| D27e-1 | 27e | `pairBy`/`pairById`, no `asList` (strict syntax, nothing left to wrap) | correction | yes | yes | `subworkflows/local/fomo_utils.nf`; lint 0/0 | valid |
| D27e-2 | 27e | optional branches gated by an empty input channel; no `when:` anywhere | improvement | yes | yes: task counts unchanged (G3), no selector warnings | `test_module_shape` (no `when:`); V1/V2 logs | valid |
| D27e-3 | 27e | MultiQC input builders not moved to a helper (retroactive, 27k) | scope cut | yes | yes | `reporting.nf`, `run_summary.nf`: two different four-element tuples | valid |
| D27e-4 | 27e | RUN_SUMMARY_TABLES (no meta) keeps tag `run_summary` (retroactive, 27k) | correction | yes | yes | `NO_META` in `tests/test_repo_conventions.py`; CLAUDE.md invariant 14 | valid |
| D27f-1 | 27f | 27b's missed A2 deletions made in 27f; correction note in 27b | correction | yes | yes | `git show --stat 61aed97` (deletions); `6b65e4b` has none | valid |
| D27f-2 | 27f | one transcript-type rule for every script | improvement | yes | yes: 23-species old-vs-new identical, G4 identical | `fomo_refactor/27f/equiv/summary.tsv` (299 rows); `tests/fomo_io` | valid |
| D27g-1 | 27g | record emitters are bash + awk (busybox-safe) | correction | yes (U6: no Python in those images) | yes | three `bin/*_record.sh` + their pytest files | valid |
| D27g-2 | 27g | the splitter writes the aggregate manifest row (`--aggregate-id`) | improvement | yes | yes | `tests/gff_by_source`; `conf/modules.config` | valid |
| D27g-3 | 27g | GFFCOMPARE_BATCH stays on its image; bash driver + awk records | correction | yes (U6: fewer image moves, none new) | yes: G4 identical; G5 identical, every `.loci` file included | the `.loci` tie-order experiment in the row; 27k B2 | valid |
| D27g-4 | 27g | batch modules loop sequentially; `process_batch` and its 123/125 retry codes deleted | improvement | yes | yes: deliverables identical; the codes came from xargs, which is gone | `modules/local/gffcompare_batch.nf`, `gff_stats.nf`; `tests/modules/test_batch_workers.py` | valid |
| D27g-5 | 27g | one `ext.identity` closure at the end of `modules.config` | improvement | yes | yes | `test_process_level_ext_follows_every_with_name` | valid |
| D27h-1 | 27h | gffstats values go through the table row builders (format, then parse back) | correction | yes (U1: `summary/tables` byte-identical) | yes | `tests/run_summary_tables` (byte-identical on V1/V2) | valid |
| D27i-1 | 27i | one GFF_STATS module; source side one task per file | improvement | yes (U5 applied) | yes: G4 identical, counts recorded | 27i Results (V2 167 → 133); stub matrix | valid |
| D27i-2 | 27i | MULTIQC has a `records` input slot; renders only with records | improvement | yes (U7) | yes | `modules/local/multiqc.nf` | valid |
| D27i-3 | 27i | render test against the 121 captured table/record pairs | correction | yes (U7) | yes | `tests/records/test_render_mqc.py`, `tests/records/fixtures` | valid |
| D27j-1 | 27j | `plans/README.md` generated by `dev/plans_index.py` | improvement | yes | yes | `tests/dev`; rerun in 27k | valid |
| D27j-2 | 27j | `plans/` check skips report text in `bin/`; the one reference removed by hand | scope cut | yes (U1: plot specs are excluded from the gate) | yes | `dev/code_metrics.py --check`; fixture diff in 27j | valid |
| D27j-3 | 27j | TD2 guard test added (27a gap) | scope addition | yes | yes: invariant 8 now has its test | `tests/modules/test_td2_predict.py` | valid |
| D27j-4 | 27j | `dev/resource_audit.py`; the subsample script takes the downloads dir | correction | yes | yes | `dev/`, `tests/dev` | valid |

No deviation needed a U-decision change, so none was escalated (master §7, step 4).

## D4. Plan coverage

Every step and every concrete sub-item of the original 27a–27j was checked against the repository
at `67ccc93` and the per-subplan commits (`git show --name-status`, the greps the plans specify,
the tests they name), not against the checkboxes.

| subplan | steps | done | replaced or cut (deviation) | superseded later | not done → fixed in 27k |
|:---|---:|---:|---:|---:|:---|
| 27a | 24 | 20 | 1 | 0 | D3 end-of-file fixer; E1 Nextflow-assets cache in CI; B1 nf-test not pinned |
| 27b | 18 | 13 | 3 | 2 | — (A2 was D27f-1) |
| 27c | 21 | 14 | 7 | 0 | — |
| 27d | 9 | 3 | 5 | 1 | — |
| 27e | 11 | 7 | 3 | 1 | B2: RUN_SUMMARY_TABLES has no meta and keeps tag `run_summary` → D27e-4 |
| 27f | 13 | 12 | 1 | 0 | B2: no direct test of `get_logger`/`script_name` |
| 27g | 12 | 5 | 5 | 2 | — |
| 27h | 7 | 6 | 1 | 0 | A3: overview, species and donors had no `emit_` step |
| 27i | 14 | 10 | 3 | 0 | C5: a dead `_mqc.tsv` publish filter on CURATE_MODELS; a comment naming `curation_to_mqc.py` |
| 27j | 12 | 8 | 4 | 0 | — |

"Superseded later" means a later subplan replaced the work on purpose: 27b A5's `transpose()` and
B1's retry codes (D27g-4), 27d A4's `process_batch` (D27g-4), 27e A5's sample-name ladders (27i),
27g D3 and E1 (27i, D27i-3). 27j C1's "3-line decision summary" is one sentence per plan
(150–250 characters), which is counted as done.

**Fixed in 27k** (each with a correction note in its subplan's Results):
- `dev/end_of_file.py` as a pre-commit hook, tested in `tests/dev`; no tracked non-data file
  needed fixing;
- an `actions/cache` step for `~/.nextflow` in the three CI jobs that run Nextflow;
- `nf-test=0.9.5` pinned in `dev/environment.yml` (the version a fresh env resolves);
- `test_get_logger_format_and_current_stderr` in `tests/fomo_io`;
- `emit_overview`, `emit_species`, `emit_donors` in `run_summary_tables.py` (the V1/V2 fixture
  test stays byte-identical);
- the CURATE_MODELS `saveAs` filter and the stale `fomo_stats.py` comment;
- D27e-4 recorded, and CLAUDE.md invariant 14 states the exception.

After these, **no step of the original plan is silently skipped.** The gates were rerun on the
clean checkout with the fixes (27k A2).

## Verdict

38 deviations (D27e-3 and D27e-4 recorded retroactively in 27k): 20 corrections (one of them a
finding), 10 improvements, 6 scope cuts (two of them moved to a later subplan and done there) and
2 scope additions. **All are valid**:
none changes a U-decision, every one keeps the invariants and the gate, and the evidence of each
is in the repository or under `/nfs/scratch01/rg/fzanarello/fomo_refactor/`. D27g-3's verdict
rests on G5, which matched B0 in all 7,666 deliverables (27k B2).

The two 27k deviations (D27k-1, D27k-2: how the plan diff and the clean checkout were made) are
procedural, change no code and are valid by the same test.

**D27k-3** (after the push, when CI ran on GitHub; it changes the gate, so it was escalated): the
golden gate compares gffcompare `.loci` files with the entries of each column sorted. Class:
correction. U1 holds: no published byte changes. The gate changed with the user's decision of
2026-10-16, recorded in 27k and master §5. Evidence:
- the same binary on the same inputs gives the golden order on NFS and another order on local
  disk or tmpfs;
- CI's files equal a local-disk rerun byte for byte;
- re-snapshotting the original pre-refactor outputs changes only the `.loci` rows (V1–V4 and
  B0);
- the refactored outputs match the new golden files on the cluster (V1–V4, G5) and on GitHub
  (V2).

Verdict: valid.
