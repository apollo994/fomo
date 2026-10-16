# 27f — Shared Python modules, filter_transcript.py, curate_models split

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: T6, T7. Review
items: D1, D6, P8, P9, N4, R5, R6, and F3/F4/F14–F19/F21/F25–F28 of the Python review.

## Status: done (2026-10-15)

## Goal

One implementation of GFF, FASTA and gzip I/O, the feature vocabulary, errors and logging,
shared by every script. The core selection step becomes tested Python instead of a two-pass awk.
No Python function exceeds complexity 15. Formatting is automatic (ruff format).

**Out of scope:** the four MultiQC adapter scripts and `run_summary_tables.py`'s input layer.
27h (run summary) and 27i (adapters) replace them; here they only switch to the shared helpers where that is trivial.

## Steps

### A. Tooling

- [x] **A1. Enable `ruff format`** (pre-commit) and format every Python file in its own step.
  The output of formatted code cannot change: G1 must pass and the golden scripts are untouched.

### B. `bin/fomo_io.py` (T7)

- [x] **B1. Write `bin/fomo_io.py`** with:
  - `open_text(path)`: gz or plain, utf-8.
  - `iter_gff(path)`: yields `(lineno, cols)` for 9-column lines, skips `#` and blank lines,
    stops at `##FASTA`.
  - `parse_attrs(col9)`: a dict. **Before choosing a duplicate-key policy, record each current
    parser's policy and pick the one that keeps every output identical.** If callers truly need
    different policies, expose a `first`/`last` parameter and say why in a docstring.
  - `iter_fasta(path)` and `fasta_names(path)`.
  - `FEATURE_VOCAB` (`GENE_TYPES`, `TX_TYPES`, `CHILD_TYPES`) with one case policy.
  - `FomoError`, `get_logger(script)`, `die()`, and `run(main)`, which maps `FomoError` to exit
    1 with one message format.
- [x] **B2. pytest for every `fomo_io` function**, including gz input, `##FASTA`, attributes
  with spaces, `Parent=a,b`, and URL-escaped values.
- [x] **B3. Prove the vocabulary is equivalent before unifying it.** (updated by D27f-2)
  - Run the old classification (curate_models' complement rule and the `TX_TYPES` set) and the
    new one on the 3 test GFF3s **and on all 23 Lycaenidae annotations** (paths in
    `samplesheet_both23.csv`).
  - Both must give the same transcript set per file. Record the comparison in Results.
  - If they differ, keep the narrower behaviour for each caller and log a deviation.

### C. Migrate the scripts to `fomo_io` (G1 + G4 after each)

- [x] **C1. `filter_gff_by_id.py`, `gff_by_source.py`:** delete the byte-identical `attr()` and
  the duplicated `TX_TYPES`; split `main()` (complexity 21 and 17) into read / select / write.
- [x] **C2. `td2_coding_filter.py`:** use `iter_fasta`, read the FASTA once, and warn (not
  silently 0.0) on a malformed PSAURON score. The output must not change.
- [x] **C3. `relocate_loci.py`:**
  - use a local `rng = random.Random(seed)` and pass it to the helpers, with the **same call
    order** so the decoys are byte-identical (G4 V2 + the pytest from 27a);
  - incremental chromosome weights **only if** they yield the identical sequence of draws;
    otherwise keep the algorithm and add a performance note to Results;
  - add a module docstring and function docstrings.
- [x] **C4. `curate_models.py` (P8):**
  - [x] split `main()` into `parse_args`, `resolve_track`, `filter_models`, `tag_ref_hits`,
    `chain_support`, `discarded_report`, `build_genes`, `novel_report`, `write_outputs`,
    `write_chains_tsv` and `write_report`;
  - [x] turn `read_reference` into a `_RefBuilder` with `add_gtf_row`/`add_gff3_row`,
    `add_fallback_intervals` and `resolve_genes`;
  - [x] use `@dataclass(slots=True)` for `Model`, `CuratedTranscript` and `CuratedGene`, a
    `ChainKey` NamedTuple, and descriptive names (no `R`/`obs`/`dis`/`nov`; no name reused for
    different things);
  - [x] replace the runtime `assert` with an explicit check; create the temp dir with
    try/finally; name the magic numbers as constants (hash length, Parent depth);
  - [x] the curation metric keys, the species buckets and `LOCATIONS` become constants in
    `fomo_stats.py`. `curation_to_mqc.py` imports them until 27i deletes it;
  - [x] the outputs stay byte-identical: the converted curate pytest suite passes, plus G4
    (curated `.gff3.gz`, `.curation.{tsv,json}`).
- [x] **C5. Unify the CLI flags of the migrated scripts (master §8).** Update the module call
  sites in the same step. List every renamed flag in Results.

### D. `filter_transcript.sh` → `filter_transcript.py` (T6)

- [x] **D1. Implement `bin/filter_transcript.py` on `fomo_io`** with the same semantics: stranded,
  spliced, `--longest` (spliced length; ties to the smallest id in byte order), funnel counts.
  The `--counts`/`--sample` row stays identical until 27g turns it into a record. It reads `.gz`
  directly.
- [x] **D2. Run the 27a `filter_transcript` tests against both implementations.** They must
  pass for both.
- [x] **D3. Real-data equivalence.**
  - Run old and new on all 23 Lycaenidae annotations × {lnc_RNA, mRNA} × {`--longest`,
    without}, as one `sbatch` job.
  - The outputs and the funnel rows must be byte-identical. Record it in Results.
  - Any difference stops the switch: fix it or log a deviation, never accept it.
- [x] **D4. Switch the `FILTER_TRANSCRIPT` module** to the Python script on `img_python`. Drop
  the scratch-decompression block and delete `bin/filter_transcript.sh` once G4 passes.

### E. Lint debt

- [x] **E1. Remove every per-file ignore from `pyproject.toml`.**
  - Fix the remaining findings: UP006/UP035 builtin generics, line length, SIM115,
    PLW1510 `check=`, unclosed handles, `encoding=`.
  - `run_summary_tables.py` and the `*_to_mqc.py` scripts may keep ignores until 27h rewrites or 27i deletes them,
    as long as each one has a comment pointing to that subplan.

## Acceptance

- G1, G2 (ruff and ruff format clean; C901 ≤ 15 everywhere except files scheduled for 27h/27i),
  G3, and G4 identical.
- The B3 and D3 equivalence reports are in Results.
- `bin/filter_transcript.sh` is gone.
- No function in the migrated scripts exceeds complexity 15; record the before and after list.

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|
| D27f-2 | B3 | if the rules differ, keep the narrower behaviour for each caller | one rule for every caller: `fomo_io.is_transcript` (anything not gene-level, structural or a child type), compared in lower case | the explicit `TX_TYPES` set and curate_models' complement rule differ only on `scRNA` (5 features in 5 of the 23 raw annotations). filter_gff_by_id and gff_by_source never read raw annotations: on all 69 files they do read in the 23-species run (46 projected, 23 filtered) both rules give the same transcripts, and their old-vs-new outputs are byte-identical (D3 harness). On a raw annotation the narrow set would orphan the children of an unlisted RNA type | none |
| D27f-1 | — (27b A2) | nothing (27b recorded these as done) | `git rm` of `assets/samplesheet_full.csv` and `assets/test_data/Lycaena_{phlaeas_282391,virgaureae_282395}/` (26 MB); `PAIRS` in `bin/subsample_test_data.sh` cut to the three test species | found while locating the 23-species samplesheet: the files were still tracked and in no samplesheet; 27b's Results claimed the deletion (correction note added there) | none (no samplesheet, test or golden file uses them) |

## Results

> **Correction (found in 27k D4):** B2 asks for a test of every `fomo_io` function; `get_logger` and
> `script_name` were covered only through `run`. 27k added
> `test_get_logger_format_and_current_stderr`.

**Equivalence on real data** (`fomo_refactor/27f/equiv.py`, one 8-CPU job, old = the scripts at
`87fd6a6`, the awk filter run in the pinned `ubuntu` image as the module ran it; summary in
`27f/equiv/summary.tsv`): every comparison byte-identical.

| check | cases | result |
|:---|:---|:---|
| filter_transcript (D3) | 23 annotations × {lnc_RNA, mRNA} × {`--longest`, none} | 92/92 identical (output and funnel row) |
| curate_models (C4) | 23 targets × {lnc_RNA, lnc_RNA.decoy}, with the samplesheet reference | 46/46 identical (`.gff3.gz`, `.merged.gff3.gz`, `.curation.{tsv,json}`); 21,019 curated genes in total |
| gff_by_source | the 46 published `allModels` GFF3s | 46/46 identical output sets |
| filter_gff_by_id, drop mode | the same 46, dropping every 7th transcript | 46/46 identical |
| filter_gff_by_id, keep mode | the 23 filtered annotations, keeping every other transcript | 23/23 identical |
| relocate_loci | 23 filtered annotations × {cap 1000, cap 0}, seed 42 | 46/46 identical (967–6,681 decoys per file) |

**B1, attribute policies before unifying.** filter_gff_by_id / gff_by_source: fields stripped,
first match; relocate_loci: fields, keys and values stripped, first match, plus a `key value`
fallback; curate_models: regex first match for ID/Parent/source/de, the reference via a dict of
stripped pairs (last match); filter_transcript.sh: no stripping, first match. A scan of the 23
annotations and 3 test GFF3s (9.95 M features) found no repeated key, no CRLF, no `##FASTA`, no
field without `=`, and whitespace only as a leading space in 121 `description=` values of gene
rows, which no script reads. So one policy (strip; first match wins; fields without `=` ignored)
keeps every output identical. **B3:** 95 files (23 raw, 3 test, 23 filtered, 46 projected): the
new rule equals curate_models' rule on 95/95 and the old `TX_TYPES` set on 90/95 (D27f-2).

**`bin/fomo_io.py`:** `open_text`, `iter_gff` (with `strict=`), `parse_attrs`, `split_ids`,
`iter_fasta`, `fasta_names`, the vocabulary (`GENE_TYPES`, `STRUCTURAL_TYPES`, `CHILD_TYPES`,
`CDS_LIKE_TYPES`, `feature_kind`, `is_transcript`), `FomoError`/`die`/`run` and `get_logger`
(one format, `<script>: <LEVEL>: <message>`, writing to the current `sys.stderr`). 10 pytest cases.

**Complexity (C901 > 10):** before — curate_models `read_reference` 38, `main` 37, `read_models`
20, `write_merged` 19, `ref_location` 12; filter_gff_by_id `main` 22; gff_by_source `main` 17;
relocate_loci `read_gff3_models` 12, `main` 12. After — curate_models `ref_location` 12 and
`write_merged` 11; nothing else in the migrated scripts above 10 (run_summary_tables
`section_aggregates` 11 is 27h's).

**C5, renamed flags (module call sites updated):** filter_gff_by_id `--keep-fasta` → `--fasta`;
td2_coding_filter `--in-fasta` → `--fasta` and `--kept-fasta`/`--coding-ids`/`--mqc` → `--prefix`
(`P.noncoding.fasta`, `P.coding_ids.txt`, `P_td2_coding_mqc.tsv`); relocate_loci `--input-gff` →
`--gff`, `--intergenic-bed` → `--bed`, `--output-gff` → `--output`; filter_transcript positional
`input TYPE > out` → `--gff --feature-type --output`. Unchanged: gff_by_source (`--outdir` +
`--name-template`: one output per source), curate_models.

**C3 performance note:** incremental chromosome weights were not added: RELOCATE_LOCI takes 1.9–2.9 s
per task in the 23-species run (cap 1000), so recomputing the weights per model costs nothing
that matters. The generator is a local `random.Random(seed)` with the old call order; a fixture
test pins the seed-7 output byte for byte.

**D:** `bin/filter_transcript.py` replaces the two-pass awk (2.6 s, 34 MB on a 440k-feature
annotation, reading `.gz` directly); FILTER_TRANSCRIPT runs it on `img_python` without the scratch
decompression. The exon type is now matched in lower case like the rest of the vocabulary
(identical on all data above).

**Other:** curation constants (`CURATION_LOCATIONS`, the table columns, the species-bucket top)
live in `fomo_stats.py`, imported by curate_models and curation_to_mqc. Two curate tests now
expect `FomoError` from a function, or `SystemExit` through `fomo_io.run`. New
`test_bin_scripts_are_executable`: the first G4 attempt failed with "Permission denied" because
the new `filter_transcript.py` lacked its execute bit (stub runs do not execute scripts).
Python lines 3,555 → 3,701 (`fomo_io.py`, docstrings and `ruff format`); pyproject keeps
per-file ignores only for the four scripts 27h/27i replace, each with its subplan.

**Gates:** G1 137 passed; G2 pass (ruff, ruff format, shellcheck, nextflow lint); G3 S1–S4 pass;
G4 V1–V4 identical. Cluster use 0.9 CPU-h.
