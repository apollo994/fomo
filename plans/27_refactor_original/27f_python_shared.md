# 27f — Shared Python modules, filter_transcript.py, curate_models split

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: T6, T7. Review
items: D1, D6, P8, P9, N4, R5, R6, and F3/F4/F14–F19/F21/F25–F28 of the Python review.

## Status: not started

## Goal

One implementation of GFF, FASTA and gzip I/O, the feature vocabulary, errors and logging,
shared by every script. The core selection step becomes tested Python instead of a two-pass awk.
No Python function exceeds complexity 15. Formatting is automatic (ruff format).

**Out of scope:** the four MultiQC adapter scripts and `run_summary_tables.py`'s input layer.
27h (run summary) and 27i (adapters) replace them; here they only switch to the shared helpers where that is trivial.

## Steps

### A. Tooling

- [ ] **A1. Enable `ruff format`** (pre-commit) and format every Python file in its own step.
  The output of formatted code cannot change: G1 must pass and the golden scripts are untouched.

### B. `bin/fomo_io.py` (T7)

- [ ] **B1. Write `bin/fomo_io.py`** with:
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
- [ ] **B2. pytest for every `fomo_io` function**, including gz input, `##FASTA`, attributes
  with spaces, `Parent=a,b`, and URL-escaped values.
- [ ] **B3. Prove the vocabulary is equivalent before unifying it.**
  - Run the old classification (curate_models' complement rule and the `TX_TYPES` set) and the
    new one on the 3 test GFF3s **and on all 23 Lycaenidae annotations** (paths in
    `samplesheet_both23.csv`).
  - Both must give the same transcript set per file. Record the comparison in Results.
  - If they differ, keep the narrower behaviour for each caller and log a deviation.

### C. Migrate the scripts to `fomo_io` (G1 + G4 after each)

- [ ] **C1. `filter_gff_by_id.py`, `gff_by_source.py`:** delete the byte-identical `attr()` and
  the duplicated `TX_TYPES`; split `main()` (complexity 21 and 17) into read / select / write.
- [ ] **C2. `td2_coding_filter.py`:** use `iter_fasta`, read the FASTA once, and warn (not
  silently 0.0) on a malformed PSAURON score. The output must not change.
- [ ] **C3. `relocate_loci.py`:**
  - use a local `rng = random.Random(seed)` and pass it to the helpers, with the **same call
    order** so the decoys are byte-identical (G4 V2 + the pytest from 27a);
  - incremental chromosome weights **only if** they yield the identical sequence of draws;
    otherwise keep the algorithm and add a performance note to Results;
  - add a module docstring and function docstrings.
- [ ] **C4. `curate_models.py` (P8):**
  - [ ] split `main()` into `parse_args`, `resolve_track`, `filter_models`, `tag_ref_hits`,
    `chain_support`, `discarded_report`, `build_genes`, `novel_report`, `write_outputs`,
    `write_chains_tsv` and `write_report`;
  - [ ] turn `read_reference` into a `_RefBuilder` with `add_gtf_row`/`add_gff3_row`,
    `add_fallback_intervals` and `resolve_genes`;
  - [ ] use `@dataclass(slots=True)` for `Model`, `CuratedTranscript` and `CuratedGene`, a
    `ChainKey` NamedTuple, and descriptive names (no `R`/`obs`/`dis`/`nov`; no name reused for
    different things);
  - [ ] replace the runtime `assert` with an explicit check; create the temp dir with
    try/finally; name the magic numbers as constants (hash length, Parent depth);
  - [ ] the curation metric keys, the species buckets and `LOCATIONS` become constants in
    `fomo_stats.py`. `curation_to_mqc.py` imports them until 27i deletes it;
  - [ ] the outputs stay byte-identical: the converted curate pytest suite passes, plus G4
    (curated `.gff3.gz`, `.curation.{tsv,json}`).
- [ ] **C5. Unify the CLI flags of the migrated scripts (master §8).** Update the module call
  sites in the same step. List every renamed flag in Results.

### D. `filter_transcript.sh` → `filter_transcript.py` (T6)

- [ ] **D1. Implement `bin/filter_transcript.py` on `fomo_io`** with the same semantics: stranded,
  spliced, `--longest` (spliced length; ties to the smallest id in byte order), funnel counts.
  The `--counts`/`--sample` row stays identical until 27g turns it into a record. It reads `.gz`
  directly.
- [ ] **D2. Run the 27a `filter_transcript` tests against both implementations.** They must
  pass for both.
- [ ] **D3. Real-data equivalence.**
  - Run old and new on all 23 Lycaenidae annotations × {lnc_RNA, mRNA} × {`--longest`,
    without}, as one `sbatch` job.
  - The outputs and the funnel rows must be byte-identical. Record it in Results.
  - Any difference stops the switch: fix it or log a deviation, never accept it.
- [ ] **D4. Switch the `FILTER_TRANSCRIPT` module** to the Python script on `img_python`. Drop
  the scratch-decompression block and delete `bin/filter_transcript.sh` once G4 passes.

### E. Lint debt

- [ ] **E1. Remove every per-file ignore from `pyproject.toml`.**
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

## Results

*(fill in)*
