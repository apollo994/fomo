# Changelog

## 0.2.0 — unreleased

A refactor: the deliverables (projected, filtered and curated annotations, alignments,
benchmarking files, run-summary tables) are byte-identical to 0.1.0 on the test profile.

### Changed
- Species ids must match `^[A-Za-z0-9_-]+$`; duplicate species abort the run at launch.
- `pipeline_info/` reports are named with the launch time (`--trace_report_suffix`), so reruns
  into the same outdir keep every report.
- MultiQC reports gain a Software Versions section (every tool, per process), list a `both`
  species' own annotation once instead of twice, and render the custom content from statistics
  records; the accuracy scatter gives every source its own colour however many there are.
- Every container image is pinned (BioContainers or Wave tags, `python` and `ubuntu` by digest)
  and declared in one place; all modules are local, trimmed copies.
- Nextflow ≥ 25.04 is required (topic channels).
- `-profile crg` honours `$NXF_WORK`.

### Removed
- The `conda` profile.

### Internal
- Statistics producers write `fomo-data/1` records; the run summary and the reports read only
  them. Per-source identity comes from a manifest, never from a file name.
- Shared Python I/O (`bin/fomo_io.py`); `bin/filter_transcript.py` replaces the awk filter.
- Tests: pytest for every script, an nf-test stub matrix, golden runs, CI; `ruff`, `shellcheck`,
  `nextflow lint` and a comment budget run in pre-commit.

## 0.1.0

The pipeline as developed through design plans 00–26: preprocessing with TD2 filtering and
optional decoys, projection with minimap2 (one alignment of all sources per target since
plans/15), all-vs-all runs with `both` species (plans/14), per-target batching (plans/18),
one filtered annotation per species with the longest isoform (plans/20), the curated annotation
(plans/21, 22, 25), no self-projection (plans/23), compressed genome I/O (plans/26), one MultiQC
report per target and a run summary (plans/16). Accuracy numbers from runs before plans/20 are
not comparable with later ones: the benchmarking reference became the longest-isoform,
TD2-filtered annotation.
