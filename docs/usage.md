# Running FOMO

The samplesheet and the roles are described in the [README](../README.md#samplesheet).

## Profiles

| profile | |
|:---|:---|
| `test` | three single-chromosome Lycaena species: one `both`, one `target`, one `source` |
| `singularity`, `docker` | the container engine |
| `crg` | the CRG cluster: SLURM queue, account, resource limits, Singularity, executor tuning |

```bash
nextflow run . -profile crg --input samplesheet.csv --outdir results
nextflow run . -profile singularity --input samplesheet.csv --outdir results --include_mrna --include_decoy
```

Run the head job inside a scheduler allocation on a cluster (it needs 1 CPU and a few GB) and
keep the work directory on a shared filesystem. With `-profile crg` the work directory is
`$NXF_WORK` when that is set.

## Parameters

`nextflow run . --help` lists them; `--help <name>` explains one. Unknown parameters are an
error, so a typo cannot silently fall back to a default.

| group | parameters |
|:---|:---|
| input/output | `--input`, `--outdir`, `--publish_dir_mode` |
| tracks | `--include_mrna`, `--include_decoy`, `--decoy_cap` (decoys per source and feature type; 0 = all), `--relocate_seed` |
| filtering | `--longest_isoform`, `--include_single_exon`, `--td2_psauron_min` |
| curation | `--curate`, `--curate_min_species` |

## Resources

Requests are set per label in `conf/base.config`; memory and time grow with each retry of a task
killed for memory or time. The genome-size-driven steps are the minimap2 index (about 10 GB per Gb
of genome) and the alignment (memory doubles per retry). Override any process for one run with a
`-c` file:

```groovy
process {
    withName: '.*:MINIMAP2_INDEX$' {
        memory = { 64.GB * task.attempt }
    }
}
```

## Rerunning

`-resume` reuses every finished task. The `crg` profile caches leniently (file size and date),
which survives the timestamp changes of a shared filesystem. Reports in `pipeline_info/` are
named with the launch time, so a rerun keeps the previous ones.

## What to look at

- `targets/<target>/curated/` — the result ([docs/output.md](output.md)).
- `targets/<target>/multiqc/multiqc_report.html` — how the target did.
- `summary/multiqc/fomo_run_summary.html` — the run across all targets and donors.
- The Software Versions section of every report lists each tool's version per process.
