# Nextflow DSL2 Template Pipeline

Minimal template for a Nextflow pipeline that follows DSL2 structure and conventions.

## Structure

- `main.nf`: entrypoint and top-level workflow
- `nextflow.config`: defaults, process settings, profiles
- `modules/local/echo/main.nf`: example process module
- `subworkflows/local/run_echo.nf`: example subworkflow using `take/main/emit`
- `assets/input/example.txt`: sample input file
- `assets/samplesheet.csv`: optional sample sheet template

## Run

```bash
nextflow run main.nf
```

Or provide your own input:

```bash
nextflow run main.nf --input 'path/to/input.txt' --outdir 'results'
```

## Customize

1. Replace `ECHO_FILE` with your real process logic.
2. Add more modules under `modules/` and compose them in `subworkflows/`.
3. Keep top-level orchestration in `main.nf`.
4. Add profile-specific settings in `nextflow.config` (e.g., Slurm, AWS Batch).
