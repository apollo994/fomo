# 27 — Codebase refactor: reduce complexity (master plan)

## Status: done (2026-10-16)

## How to run it

Trigger the whole implementation with this single prompt:

```
Implement plans/27_refactor.md — follow its "Execution protocol" from the first subplan that is not "done".
```

The same prompt resumes after an interruption or a context reset, because the progress lives in
this file and in the subplans. Do not put extra instructions in the prompt: if something needs to
change, change the plan first.

---

## 1. Why

The code review ([`ablation/CODE_REVIEW.md`](../ablation/CODE_REVIEW.md), 2026-10-15) found three
root causes behind most of the maintenance cost.

1. **Identity is carried in strings and recovered by parsing.**
   - The track name is assembled inline 35 times.
   - Per-source filenames are parsed by three parsers.
   - MultiQC display labels (`1_raw`, `2_lncRNA`) are parsed back as data.
2. **Prose stands in for structure and tests.**
   - The Nextflow layer has about one comment line per code line (0.85).
   - CLAUDE.md is 514 lines.
   - At least 20 comments contradict the code.
   - Only one script has a test, there is no CI, and the stub run does not complete.
3. **nf-core form without nf-core tooling.**
   - Two "nf-core" modules are forks, and two are unused.
   - Software versions are emitted but never collected.
   - Resource labels exist only to be overridden.
   - The conda profile works for only 3 of 21 local modules.

The **goal is lower complexity at identical results**: fewer concepts, fewer copies, fewer
processes and fewer lines. Every invariant becomes enforced by a test instead of explained in
prose. Finding IDs (P1, D3, N1, C1, V1, R1, T2, …) refer to `ablation/CODE_REVIEW.md`.

## 2. Goals and non-goals

**Goals** (measured in 27k):

| metric | before | target | after (27k, measured) |
|:---|---:|---:|---:|
| comment : code ratio, Nextflow workflows + config | 0.83 / 0.88 | ≤ 0.40 | 0.26 / 0.20 |
| run / incident / history phrases in code (`code_metrics.py`) | 78 | 0 | 0 (budget lists); 5 with the review's regexes, see below |
| `CLAUDE.md` | 514 lines | ≤ 220 lines | 150 lines |
| copies of the track-name expression | 35 | 0 (`meta.track`) | 0 |
| parsers of the per-source filename grammar | 3 + 1 partial | 0 (manifest) | 0 |
| MultiQC adapter processes (`*_TO_MQC`, `GFFCOMPARE_TO_SCATTER`) | 5 | 0 | 0 |
| GFF statistics modules | 3 | 1 | 1 (aliases `GFF_STATS_TARGET`, `GFF_STATS_PROJECTED`) |
| `modules/nf-core` | 9 modules, 2 forked | gone | gone |
| tests | 1 script, not pytest | pytest for every bin script + nf-test stub matrix + golden runs + CI | 277 pytest, nf-test S1–S4, golden V1–V4 + real-data snapshot, CI (4 jobs) |
| `nextflow lint` errors | 6 | 0 | 0 (and 0 warnings, from 15) |
| Python functions above complexity 15 | 6 | 0 | 0 (maximum 12, from 38) |
| tasks in variant V2 (test, mRNA + decoy) | baseline from 27a | measurably fewer (recorded) | 167 → 133 (V1 61 → 47, V3 59 → 46; real data 1,154 → 878) |

Every goal is met. The phrase row: the budget check (`dev/code_metrics.py --check`) uses 27b's
phrase lists and finds nothing; the review's broader regexes still match five present-tense
lines, four saying what a record "measured" and one docstring saying the gffcompare cliques are
"verified against" the chains derived in Python (D27b-4).

**Non-goals:**

- **No scientific behaviour changes.** The curation rule, donor selection, TD2 thresholds,
  `splice:hq` and the mRNA coding mask recommended in `ablation/ABLATION_REPORT.md` are future
  plans.
- No new container images (U6).
- No change to the published directory layout or deliverable file names. The only exception is
  `pipeline_info/` (U5).
- No push. No change to remote state.

## 3. Decisions

**User decisions (2026-10-15).** These are fixed. A deviation that touches one of them is a
**stop condition** (§7).

| id | decision |
|:---|:---|
| **U1** | **Regression gate: deliverables identical.** Published GFF3/BAM/gffcompare/curated/filtered files and `summary/tables` values must be identical on the golden variants (§5). MultiQC HTML, MultiQC data, section layout and `pipeline_info/` may change. |
| **U2** | **All modules local; no `modules/nf-core` directory and no `modules.json`.** Modules are trimmed to what FOMO calls. nf-core remains a source of prior art, not a dependency. |
| **U3** | **Tests: pytest + nf-test + GitHub Actions (docker) + pre-commit**, built first, as the safety net (27a). |
| **U4** | **Software versions: collect the `versions` topic into a MultiQC section.** Fill the modules that report nothing and fix the broken evals. |
| **U5** | **Allowed behaviour changes, and only these:** <br>• species pattern `^[A-Za-z0-9_-]+$`; <br>• timestamped `pipeline_info` files; <br>• no duplicate target-side GFF stats for `both` species. <br>**`--include_single_exon` stays.** |
| **U6** | **Containers: pin the existing images (digest or immutable tag), declared in one config file, and drop the conda profile.** No new images. |
| **U7** | **Stats architecture: producers emit data records (JSON).** MultiQC content is rendered from records, and the run summary reads only records. |
| **U8** | **Commits:** commit `ablation/` on branch `ablation`, then create branch `refactor`, then one commit per finished subplan whose message references the plan. **Never push.** |

**Technical decisions.** These were taken by the plan and may be revised through the deviation
protocol (§7) as long as no U-decision is affected.

| id | decision | why |
|:---|:---|:---|
| T1 | **Data record format `fomo-data/1`**: `{"schema":"fomo-data/1","section":…,"identity":{…},"data":{…}}`, one file `<prefix>.<section>.data.json` per item (a list of records is allowed). Identity keys are drawn from `species, side, target, source, track, label, order, kind` (and `stage`, as 27g specifies). Payloads are tool-specific; for gff-feature-stats the payload is its JSON verbatim. | U7 |
| T2 | **The MultiQC renderer (`bin/render_mqc.py`) runs inside the per-target MULTIQC task.** That image has Python 3.14. The tool-only images (samtools, seqkit, gffread, bedtools, minimap2) have no Python, so samtools and seqkit records are emitted by two small tested awk scripts in `bin/` (bash wrappers around awk, D27g-1). `GFFCOMPARE_BATCH` moves to the CURATE image, which has Python and gffcompare 0.12.6, the same version. (updated by D27g-3: it stays on its BioContainers image, with a bash driver and an awk record emitter, because gffcompare's `.loci` tie order depends on the build.) | U6 forbids new images; checked 2026-10-15 by `singularity exec` on the cached images |
| T3 | **Per-source identity comes from a manifest** written by the splitter (`gff_by_source.py` → `manifest.tsv`). Filenames are never parsed. | P1 |
| T4 | **Track registry in one Groovy function file**, `subworkflows/local/fomo_utils.nf`. Meta keys `track` (dot form, e.g. `lnc_RNA.decoy`), `label`, `order` and `td2` are stamped once where tracks are created. Config closures read meta keys only. | N1, P2 |
| T5 | **Config layout:** nf-core include order (base → profiles → modules); selectors `'.*:NAME$'`; resources only in `conf/base.config`, with labels that state real needs; `conf/crg.config` holds site settings only; containers in `conf/containers.config`, keyed by one `img_*` label per module. | P5, P6, C4, C5 |
| T6 | **`filter_transcript.sh` is rewritten as `filter_transcript.py`** on the shared I/O module, with byte-identical output. | P9 |
| T7 | **Shared Python modules as flat siblings in `bin/`:** `fomo_io.py` (I/O, GFF/FASTA parsing, errors, logging) and `fomo_stats.py` (naming, identity, record helpers, gffcompare parser). They are imported with the existing `sys.path` pattern. A package directory was rejected because it is unverified across executors. | D1, D2 |
| T8 | **Comment policy** (§8) applies to every line written or touched. | V1–V5 |
| T9 | **Golden regression** with `tests/golden/` (snapshot + compare) on four test-profile variants. Real-data check (G5) in 27k. | U1 |
| T10 | **Real runs on this cluster:** `sbatch` a job (4 cpu, 16 GB) that runs `-profile test,singularity` with the local executor inside it. CI uses `-profile test,docker`. Scratch root: `/nfs/scratch01/rg/fzanarello/fomo_refactor/`. | docker is not available on the cluster; the interactive job has 1 cpu |

### Precedence over CLAUDE.md during the refactor

CLAUDE.md is loaded as project instructions in every session, and it is rewritten only in 27j.
Until then, where it conflicts with this plan, **this plan wins**: it was approved after
CLAUDE.md, with decisions U1–U8. Known conflicts:

- "prefer nf-core module conventions / modules" — U2 replaces it with local modules;
- "never nest subworkflows", which exists only for the fully-qualified selectors that 27d
  removes;
- the `ext.when` gating of curation (27e-C2);
- the descriptions of the adapter modules and `*_mqc.tsv` routing (27g–27i);
- the gunzip, filter and other wording that 27b–27f change.

The **data contracts** in CLAUDE.md are not conflicts: they are the invariants in §4 and must
hold throughout. Keep the Seqera-MCP guidance: search nf-core for prior art, then vendor the
module.

## 4. Invariants that must survive

These are the contracts in today's CLAUDE.md. Each must be enforced by a test by the end of 27k;
the subplan that owns the test is in brackets.

1. **Filtered-annotation contract.** `<source>.<ft>.filtered.gff3` describes exactly the
   projected FASTA, matched by exact id, and keep-mode fails on mismatch. [27a pytest, 27f]
2. **Merge/split custody.** `>tid|type|species` → BAM QNAME → `source=` → the split. Projected
   IDs are `<tid>|<source>`, and a duplicate id fails loudly. [27a, 27b]
3. **No self-projection.** A `both` target's records are dropped before alignment, on every
   track. [27a stub/golden, 27c]
4. **No gunzip task and no faidx shortcut.** Inputs stay compressed; FASTA input is plain gzip in
   production. [27a variant V4]
5. **Decoy placement is annotation-wide** (raw GFF3 → gene BED) and reproducible with the seed.
   [27a golden V2]
6. **Feature types are derived once.** The projection ⋈ reference join keys on
   `[target, feature_type]`. [27a stub matrix]
7. **A target with nothing projected yields an honest empty result** (`remainder: true` path).
   [27a stub or pytest]
8. **TD2 guards:** an empty FASTA, and zero complete ORFs, mean "nothing filtered". [27a pytest]
   (updated by D27j-3: 27a did not write it; `tests/modules/test_td2_predict.py` was added in 27j.)
9. **`FILTER_TRANSCRIPT` rules:** stranded, spliced, longest isoform, ties broken by the
   smallest id, and the funnel counts. [27a pytest → oracle for 27f]
10. **Curated IDs are content-hashed and stable.** [existing curate tests → pytest in 27a]
11. **One MultiQC report per target plus one aggregate-only run summary**, and the publish
    layout. [27a stub matrix (published paths)]
12. **The task-count formulas** (F·S, F·D·T, …). [27a stub matrix]

## 5. Regression gate

Every subplan ends by running the gates listed in its **Acceptance** section.

| gate | command (from 27a on) | when |
|:---|:---|:---|
| **G1 unit** | `pytest -q` in the dev environment (`dev/environment.yml`) | every subplan, often |
| **G2 lint** | `pre-commit run --all-files`: ruff, shellcheck, `nextflow lint` (0 errors), and from 27j the `code_metrics` budget | every subplan |
| **G3 stub** | `nf-test test tests/pipeline --profile +singularity` (stub matrix S1–S3) | every subplan |
| **G4 golden** | `tests/golden/run_variant.sh V1..V4` (via `sbatch`, T10), then `tests/golden/compare.py`. **Zero differences in deliverables**, and V4 must equal V1. | every subplan |
| **G5 real data** | 23 Lycaenidae `both` species, `--include_decoy`, `-profile crg`, snapshot compared with `tests/golden/realdata_B0.tsv` | 27k (optionally earlier) |

**Golden variants** (all `-profile test`):

| variant | definition |
|:---|:---|
| V1 | default (lncRNA only) |
| V2 | `--include_mrna --include_decoy` |
| V3 | `assets/samplesheet_target_noref.csv` |
| V4 | V1 with plain-gzip FASTAs and uncompressed GFF3s, generated on the fly; must equal V1 |

**What "identical deliverables" means** (normalisation implemented in `tests/golden/snapshot.py`):

| published path | compared as |
|:---|:---|
| `targets/*/annotation/*.gff3`, `sources/*/annotation/*.gff3` | bytes |
| `targets/*/alignment/*.bam` | `samtools view` records (no header) |
| `targets/*/gffcompare/*.stats` | lines not starting with `#` (the header holds the command line) |
| `targets/*/gffcompare/*.loci` | the entries of each column sorted: gffcompare orders tied entries by heap layout, which follows the filesystem's block size (D27k-3, user decision 2026-10-16) |
| other `targets/*/gffcompare/*` | bytes |
| `targets/*/curated/*.gff3.gz` | decompressed bytes |
| `targets/*/curated/*.curation.tsv` | bytes |
| `targets/*/curated/*.curation.json` | canonical JSON, without `params.gff` / `params.ref` (staging paths) |
| `summary/tables/run_summary.json` | canonical JSON |
| `summary/tables/run_*_mqc.tsv` | bytes |
| `*.csi`, `summary/tables/*_mqc.json` (plot specs), `*/multiqc/**`, `pipeline_info/**` | excluded |
| **the set of paths itself** | identical |

A failing G4 is never fixed by updating the golden files. The only exception is a difference
that a U-decision explicitly allows, and none of U5's changes touches a deliverable.

## 6. Subplans

| id | file | title | depends on | review items | est. | status |
|:---|:---|:---|:---|:---|---:|:---|
| 27a | [27a_safety_net.md](27a_safety_net.md) | Golden baseline, pytest, nf-test, lint, CI | — | T1–T3, R4 | 3 d | done (2026-10-15) |
| 27b | [27b_delete_and_fix.md](27b_delete_and_fix.md) | Delete dead code, fix latent bugs, comment sweep | 27a | §2, R1–R3, V1–V4 | 2 d | done (2026-10-15) |
| 27c | [27c_local_modules.md](27c_local_modules.md) | All modules local, containers, versions, no conda | 27b | C1–C5, D4 (filters) | 3 d | done (2026-10-15) |
| 27d | [27d_config_structure.md](27d_config_structure.md) | Config layout, selectors, resources, pipeline_info | 27c | P5, P6, R7, C4 | 1.5 d | done (2026-10-15) |
| 27e | [27e_track_identity.md](27e_track_identity.md) | Track registry, meta identity, channel simplification | 27d | N1, N2, P7, P10 | 2.5 d | done (2026-10-15) |
| 27f | [27f_python_shared.md](27f_python_shared.md) | Shared Python modules, filter_transcript.py, curate split | 27e | D1, P8, P9, N4, R5, R6 | 4 d | done (2026-10-15) |
| 27g | [27g_stats_records.md](27g_stats_records.md) | Producers emit data records (additive); manifest | 27f | U7, T1–T3, P1, D3 | 2.5 d | done (2026-10-15) |
| 27h | [27h_run_summary_records.md](27h_run_summary_records.md) | Run summary reads records | 27g | P2, P4, R5 | 2 d | done (2026-10-15) |
| 27i | [27i_multiqc_from_records.md](27i_multiqc_from_records.md) | Renderer + section registry; delete adapters, one GFF_STATS | 27h | P3, D2, D4, U5 | 3 d | done (2026-10-15) |
| 27j | [27j_docs.md](27j_docs.md) | CLAUDE.md, README, docs/, ADR plans, CHANGELOG, budget | 27i | V5–V7, T4, T5 | 2 d | done (2026-10-15) |
| 27k | [27k_final_validation.md](27k_final_validation.md) | Clean-checkout gates, real-data G5, plan comparison | all | — | 1.5 d | done (2026-10-16) |

The estimates are engineering days and are only for ordering. The order matters: 27a protects
everything after it. 27c–27e restructure without changing logic. 27f–27i change the code paths
under the golden gate. Records are added (27g) before any consumer switches to them
(27h run summary, 27i MultiQC), and consumers switch before the old path is deleted (27i).

## 7. Execution protocol

### Step 0 — once, before 27a

1. Check: the current branch is `ablation`, and the working tree is clean except for
   `ablation/` and `plans/27*`.
2. Commit `ablation/` on `ablation` with the message `ablation study: scientific and code
   reviews`.
3. Run `git switch -c refactor`.
4. **Freeze the original plan:** if `plans/27_refactor_original/` does not exist, copy every
   `plans/27*.md` into it. **Never edit those copies afterwards.**
5. Commit `plans/27*` and the frozen copy with the message `plans: 27 refactor master plan
   (+ frozen original)`.
6. Set this file's status to `in progress`.

### For each subplan, in table order, starting at the first one not marked `done`

1. **Read** this master plan (§3–§8) and the subplan in full. Set the subplan status to
   `in progress` here and in the subplan.
2. **Implement** the steps in order and tick each checkbox as it is done. Run G1 and G2 on
   touched files often. Long runs go through `sbatch` (T10); poll them with a background
   until-loop, not with sleeps.
3. **When reality differs from the plan, apply the deviation protocol below before continuing.**
4. **Accept:** run every gate in the subplan's Acceptance section. Record the commands and their
   results in the subplan's **Results** section, with numbers.
5. **Update:** the subplan status becomes `done (YYYY-MM-DD)`, the table in §6 is updated, and
   the **Status log** (§10) gets a line.
6. **Commit** once:

   ```
   refactor(27x): <one-line summary>

   Plan: plans/27x_<name>.md (master: plans/27_refactor.md)
   Deviations: D27x-1, …  (or "none")
   ```

   Use this session's commit attribution trailer. **Never push.**
7. **Continue** with the next subplan without waiting, unless a stop condition holds.

### Deviation protocol

Use it when a step is wrong, impossible, insufficient, or a better route is found.

1. **Fix the working plan.** Edit the subplan steps so they describe what is actually done; the
   working plan always matches the code.
2. **Log it.** Append a row to the subplan's **Deviations** table:

   ```
   | D27x-N | step | original plan said | what was done | why (evidence: command, output, file:line) | impact on later subplans |
   ```

3. **Propagate.** If later subplans are affected, edit them and add a note `(updated by D27x-N)`
   at the changed step.
4. **Escalate.** If the deviation would change a U-decision, the gate (§5) or an invariant (§4),
   **stop and ask the user** before doing it.
5. **Never edit `plans/27_refactor_original/`.**

### Stop conditions

Write the current state into the subplan, commit the work in progress *only if every gate
passes*, then ask the user, when any of these holds:

- (a) The work conflicts with a U-decision or an invariant.
- (b) A gate fails and the cause is not understood after a reasonable investigation.
- (c) G4 or G5 shows a deliverable difference that is not a proven bug fix. **Never regenerate
  golden files to make a gate pass.**
- (d) Cluster use would exceed about 40 CPU-hours in total for the refactor. Track it in §10.
- (e) Credentials, a push, or any action outside this repository and the scratch root would be
  needed. Installing tools into a dev conda environment under `$HOME` is allowed.

### Resuming

After an interruption or context reset, the §6 table, §10 and the subplans' checkboxes are the
source of truth. Re-read this file and continue at the first unchecked step of the first
subplan that is not `done`. Check `git status` and `git log` first: an interrupted step may have
left changes that are not committed.

### Final comparison

27k compares this working plan with `plans/27_refactor_original/` and validates every deviation.
The result goes to `plans/27z_deviation_review.md`.

## 8. Conventions for all subplans

**Comments (T8).**
- Present tense: they say what is true now and why. At most about 3 lines, unless they state an
  invariant.
- **No** run names, job IDs, clades, measurements or incident stories.
- **No** history ("used to", "no longer", "old X", "since plans/NN").
- **No** `file:line` references, and no names of things that no longer exist.
- A `plans/NN` pointer is allowed only as `Rationale: plans/NN`, only when that plan is the
  current rationale, and never in user-facing text (schema help, MultiQC text, README).
- History goes in commit messages, CHANGELOG and plans.

**Naming.**
- Processes are `VERB_OBJECT`; a module file name is the process name in lower case.
- Emits: `mqc` (MultiQC files), `data` (records). Versions use the topic only.
- Feature types use Sequence Ontology terms in identifiers and filenames (`lnc_RNA`, `mRNA`);
  "lncRNA" appears in prose and labels only.
- Every track is `meta.track` (T4).

**Python CLI.**
- Inputs `--gff/--fasta/--json`; outputs `--output` (one file) or `--prefix` (several).
- `--sample` for the row name; booleans are `store_true`.
- Errors raise `FomoError`; `main()` turns them into exit 1 with one logger format.

**Shell.** Every script starts with `set -euo pipefail` and `export LC_ALL=C`, and must be
shellcheck-clean.

**Modules.**
- One `img_*` label, no `container` or `conda` lines (U6, T5).
- Tag: `"${[meta.target_id, meta.id, meta.track].findAll().unique().join('.')}"` or the
  subplan's refined form.
- A stub is mandatory. `when:` only where config gates the module.

**Tests.** Every new or changed script gets pytest coverage in the same subplan. A contract shared
by two files (a filename, a label, a section id) gets a consistency test.

## 9. Risks

| risk | mitigation |
|:---|:---|
| Golden runs are not deterministic (threads, PSAURON, sort ties) | 27a runs every variant twice before anything changes; normalise only proven non-determinism, and document it in 27a |
| CI cannot be verified because nothing is pushed | validate the workflow YAML and run every CI command locally (Singularity instead of Docker); state it in 27j's Outcome |
| Docker images behave differently from Singularity in CI | same image digests (U6); the golden comparison in CI is the check |
| Local tool versions differ from the pipeline's (local gffcompare 0.12.10 vs pinned 0.12.6) | `dev/environment.yml` pins the pipeline versions for integration tests |
| `filter_transcript.py` differs from the awk on real annotations | 27f compares old against new on all 23 Lycaenidae annotations, not only test data |
| Long autonomous run loses context | the protocol is resumable (§7); commit per subplan; small checkboxes |
| Cluster cost | budget in stop condition (d); G5 only once (27k) |

## 10. Status log

| date | subplan | event | CPU-h used (cumulative) |
|:---|:---|:---|---:|
| 2026-10-15 | — | plan written; decisions U1–U8 recorded | 0 |
| 2026-10-15 | 27a | done: golden baseline (deterministic, V4 = V1), 81 pytest, stub matrix S1–S3, lint clean, CI written; 4 deviations | 1.5 |
| 2026-10-15 | 27b | done: −1,455/+496 lines, 4 latent bugs fixed, comment ratio of the Nextflow layer 0.85 → 0.21, G4 identical; 4 deviations | 2.1 |
| 2026-10-15 | 27c | done: `modules/nf-core` and conda gone, 27 local modules (1,272 → 964 lines), 12 `img_*` labels in `conf/containers.config` (python/ubuntu pinned by digest), versions topic → Software Versions in every MultiQC report, G4 identical; 7 deviations (D27c-6: floor → 25.04 in 27d) | 3.0 |
| 2026-10-15 | 27d | done: include order base → containers → profiles → modules, crg site-only, resources only in `base.config` (crg values identical for all 31 processes), 38 `.*:NAME$` selectors, TD2 told apart by `meta.stage`, timestamped `pipeline_info`, floor 25.04, G4 identical; 4 deviations | 3.5 |
| 2026-10-15 | 27e | done: track keys stamped once (`withTrack`), no pseudo meta ids, one tag rule, 9 pairings → `pairBy`, optional branches gated by empty input (no `when:`, no selector warnings), lint warnings fail the hook, G4 identical; 2 deviations | 3.9 |
| 2026-10-15 | 27f | done: `fomo_io.py` shared by every migrated script, curate_models split (max complexity 38 → 12), `filter_transcript.py` replaces the awk, flags unified; 299/299 old-vs-new comparisons on the 23-species data byte-identical, G4 identical; 27b's missed A2 deletions made; 2 deviations | 4.8 |
| 2026-10-15 | 27g | done: every stats producer also writes fomo-data/1 records (7 sections), a manifest carries per-file identity, the bridge test proves records hold every table value (121 pairs), G4 identical; GFFCOMPARE_BATCH stays on its image (gffcompare `.loci` tie order depends on build and path); 5 deviations | 6.7 |
| 2026-10-15 | 27h | done: the run summary reads only data records (cached indexes, pure compute + emit, no label or filename parsing, no lint ignores), byte-identical on the V1/V2 fixtures and in G4; 1 deviation | 7.0 |
| 2026-10-15 | 27i | done: MultiQC custom content rendered in MULTIQC from records via one registry; 5 adapter processes, 4 scripts and the producers' tables deleted; one GFF_STATS module; duplicate target rows gone (U5); V2 167 → 133 tasks; G4 identical; 3 deviations | 7.3 |
| 2026-10-15 | 27j | done: CLAUDE.md 514 → 150 lines with 16 tested invariants, README + docs/, CHANGELOG, CONTRIBUTING, decision-record headers on plans 00–26, dev/ tools with the comment budget and link check in pre-commit; TD2 guard test added; G4 identical; 4 deviations | 7.6 |
| 2026-10-16 | 27k | done: G1–G4 pass on a clean checkout with a fresh dev env; G5 (23 species, real data) identical to B0 in all 7,666 deliverables with 878 tasks instead of 1,154; goals measured; 65 plan hunks classified, 38 deviations valid (2 retroactive), 8 unfinished items found by the coverage audit and finished; 2 deviations | 15.6 |

## 11. Outcome

**Result.** The refactor reached every goal in §2 with identical deliverables. Every gate passes
on a clean checkout with a dev environment built fresh from `dev/environment.yml`: 277 pytest,
lint, the nf-test stub matrix S1–S4, and golden V1–V4 identical. On real data (G5: the 23
Lycaenidae `both` species, lncRNA + decoy, `-profile crg`), all 7,666 deliverables are identical
to the pre-refactor run B0: 0 missing, 0 extra, 0 changed, `.loci` files included.

**Goals.** The table in §2 has the measured values; every target is met.

**Size.**
- Pipeline code (`main.nf`, `workflows`, `subworkflows`, `modules`, `conf`, `nextflow.config`,
  `bin`, `assets/multiqc`): 7,740 → 6,087 lines (`git diff c86521c`: +4,102 / −5,743).
- Modules: 30 files, 1,563 lines → 22 files, 763 lines; 30 → 22 processes.
- Code and comment lines per layer, before → after:

  | layer | code | comment |
  |:---|---:|---:|
  | workflows + subworkflows | 452 → 358 | 373 → 94 |
  | modules | 1,189 → 621 | 243 → 43 |
  | config | 415 → 351 | 366 → 71 |
  | `bin/` Python | 2,467 → 2,630 | 492 → 376 |
  | `bin/` shell | 341 → 273 | 102 → 61 |
  | MultiQC YAML | 355 → 352 | 54 → 40 |
  | tests (code only) | 467 → 2,075 | 44 → 124 |
  | `dev/` | — → 413 | — → 55 |

  `bin/` Python grows because of `fomo_io.py`, the records, the renderer and docstrings, while
  the adapter scripts and five processes are gone. Test fixtures add about 37,000 data lines.
- Lint: ruff over `bin`, `dev` and `tests` with the 27a rule set went from 251 findings to 0, with
  no per-file ignores; `nextflow lint` from 6 errors and 15 warnings to 0 and 0.

**Tasks.** V1 61 → 47, V2 167 → 133, V3 59 → 46. On the 23-species run, 1,154 → 878 (−24 %):
the five adapter processes and the duplicate target-side statistics (U5) are gone, and every
other process runs the same number of tasks.

**Cluster use.** 15.6 CPU-h in total (27k: 8.0, most of it G5), within the 40 CPU-h budget.

**Deviations** ([27z_deviation_review.md](27z_deviation_review.md)). There are 38 in 27a–27j,
plus 2 procedural ones in 27k. By class: 20 corrections, 10 improvements, 6 scope cuts and 2
scope additions. All are valid, and none touched a U-decision.

Findings worth remembering:
- **27b A2 slip.** 27b recorded deletions it had not made. They were found and made in 27f
  (D27f-1).
- **TD2 guard test.** The test 27a owed for invariant 8 was added in 27j (D27j-3).
- **gffcompare `.loci` tie order** depends on the gffcompare build and the run directory
  (D27g-3), so GFFCOMPARE_BATCH stays on its image. G5 shows no such difference.
- **Coverage audit (27k).** It found eight small unfinished items in 27a, 27e, 27f, 27h and 27i,
  all finished in 27k:
  - the end-of-file hook;
  - the CI Nextflow cache;
  - the nf-test pin;
  - the `get_logger` test;
  - three missing run-summary emit steps;
  - a dead publish filter and a stale comment;
  - the RUN_SUMMARY_TABLES tag exception (D27e-4).
- **Propagation misses.** Four places in the plan text did not yet show a recorded deviation.
  They were fixed in 27k.

**Not verified, because nothing was pushed:**
- the GitHub-hosted runners;
- `setup-micromamba` and its cache;
- `setup-nextflow`;
- the `~/.nextflow` cache step;
- artifact upload;
- Docker execution of the pinned images, including the `linux/amd64` pin on arm64 hosts.

Every CI `run:` command was executed locally with Singularity in place of Docker.

**Follow-ups (candidate plans 28+).** These are the scientific recommendations of
`ablation/ABLATION_REPORT.md` §5, which this refactor left out by design (§2 non-goals):
1. **An adaptive donor set:** rank donors by median `de` and curate over the nearest N, with all
   donors as the fallback. Also publish the nearest donor's projection as a high-recall
   companion.
2. **Per-track curation thresholds** (`curate_min_species` per track).
3. **An mRNA coding mask** for unannotated targets, which would make `--include_mrna` part of
   their deliverable.
4. **A `splice:hq` minimap2 option,** or an automatic choice for close donors.
5. **A decoy FDR estimator** that sends the decoy and real tracks through the same reference
   filter, or the support calibration in its place.
6. Smaller items:
   - re-document or remove `--include_single_exon`, which U5 kept;
   - review the `antisense_exonic` class.

**Branch.** `refactor` holds one commit per subplan (27a–27k). It is not pushed and is ready for
review.
