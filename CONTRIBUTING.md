# Contributing

## Principles

- **One home per fact.** Behaviour and invariants live in code and tests; a comment says why in a
  few lines; [CLAUDE.md](CLAUDE.md) indexes the invariants and conventions; `docs/` explains the
  pipeline to people; `plans/` records decisions; the [CHANGELOG](CHANGELOG.md) records history.
- **Contracts are tested.** A rule that two files must agree on (a file name, a label, a section
  id) has a test that fails when they drift.
- **Deliverables do not move silently.** A change to a published file is a decision, made in a
  plan and visible in the golden gate.
- **Small, plain code.** No abstraction without a second user; delete what nothing calls.

## Conventions

The conventions (comments, naming, modules, configuration, Python, shell) and the checklists for
adding a track, a statistics section, a parameter or a module are in [CLAUDE.md](CLAUDE.md).

## Gates

Every change passes, in order:

1. `pytest -q`
2. `pre-commit run --all-files` — ruff, ruff format, shellcheck, `dev/nextflow_lint.sh` (errors
   and warnings), `dev/code_metrics.py --check` (comment budget), `dev/check_links.py`, trailing whitespace and
   `dev/end_of_file.py`
3. `nf-test test tests/pipeline --profile +singularity` — the stub matrix (task counts,
   published paths)
4. `tests/golden/gate.sh <scratch dir>` — the golden variants V1–V4; deliverables must be
   identical. A difference is investigated, never accepted by regenerating golden files.

The development environment is `dev/environment.yml` (`conda env create -f dev/environment.yml`).

## Running CI locally

CI (`.github/workflows/ci.yml`) runs the same four gates with Docker. Locally, with the
development environment active, run the gates above; the golden gate needs SLURM (or
`tests/golden/run_variant.sh V1 <dir> test,docker` per variant, then `tests/golden/compare.py`).

## Commits

One commit per finished change, with a subject that says what changed (`area: change`), a body
that says why, the plan it implements (`Plan: plans/NN_...`) and any deviation from it.
