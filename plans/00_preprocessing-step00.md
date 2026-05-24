# Plan: FOMO Preprocessing — Step 00 (incremental)

## Context

FOMO transfers lncRNA annotations from source species onto a target assembly. The preprocessing stage prepares spliced FASTAs that the projection stage (minimap2) consumes. This plan covers **only step 00** ("extract lncRNA/mRNA spliced FASTA") plus the project skeleton needed to run it end-to-end on subsampled test data. Steps 02–05 (intergenic intervals, decoy relocation, decoy FASTA, statistics) come in follow-up work.

## Decisions

| Decision | Choice |
|----------|--------|
| Input model | Single samplesheet CSV with `role` column (`target` / `source`) |
| Scope | Step 00 only (incremental) |
| Container strategy | nf-core modules where available; prefer gffread over AGAT |
| Feature types | Both `lnc_RNA` and `mRNA` |
| Longest isoform filtering | Skipped (dropped per gffread-only decision) |
| Project layout | nf-core style |
| Test target species | `Lycaena_phlaeas_282391` |
| Subsampled test data location | Replace symlinks in `assets/test_data/` with committed files |
| Subsampling script location | `bin/subsample_test_data.sh` |
| Module shape | Two modules: `FILTER_TRANSCRIPT` (local) + `GFFREAD` (nf-core) |

## 1. One-off: Subsample test data to one chromosome per species

**Script:** `bin/subsample_test_data.sh` — run once by hand. Replaces symlinks under `assets/test_data/<species>/<assembly>/` with committed subsampled files.

Per species dir:
1. `gunzip -c <species>.fna.gz | samtools faidx -` → write `.fai` to scratch.
2. Read `.fai`, pick the contig with the largest value in the length column → `LONGEST`.
3. `samtools faidx <species>.fna.gz "$LONGEST" | bgzip > <species>.fna.gz.new`
4. `zcat <species>.gff3.gz | awk -v c="$LONGEST" -F'\t' '/^#/ || $1==c' | gzip > <species>.gff3.gz.new`
5. Delete the symlinks; move `.new` files into place.

After running: update `.gitignore` to remove the `assets/` exclusion so these files are committed.

> Tools: `samtools`, `bgzip`, `awk`, `gzip`. Dev utility only — not containerized.

## 2. Project skeleton (nf-core layout)

New files to create:

```
main.nf
nextflow.config
nextflow_schema.json
workflows/fomo.nf
subworkflows/local/preprocessing.nf
modules/local/filter_transcript.nf
modules/nf-core/gffread/main.nf           # installed via `nf-core modules install gffread`
bin/filter_transcript.sh                  # imported verbatim from /Users/fzanarello/work/repo/dirty_scripts/annotation/filter_transcript.sh
bin/subsample_test_data.sh
conf/test.config
conf/modules.config
assets/samplesheet.csv
```

## 3. Samplesheet (`assets/samplesheet.csv`)

```csv
species,role,fasta,gff3
Lycaena_phlaeas_282391,target,assets/test_data/Lycaena_phlaeas_282391/GCA_905333005.2/Lycaena_phlaeas_282391.fna.gz,assets/test_data/Lycaena_phlaeas_282391/GCA_905333005.2/Lycaena_phlaeas_282391.gff3.gz
Lycaena_alciphron_282377,source,assets/test_data/Lycaena_alciphron_282377/GCA_964273875.1/Lycaena_alciphron_282377.fna.gz,assets/test_data/Lycaena_alciphron_282377/GCA_964273875.1/Lycaena_alciphron_282377.gff3.gz
Lycaena_hippothoe_580924,source,...,...
Lycaena_thersamon_265360,source,...,...
Lycaena_virgaureae_282395,source,...,...
```

In `workflows/fomo.nf`:
- Parse with `splitCsv(header:true)`, branch into `target` and `source` channels by `role`.
- Step 00 only needs `source` rows — target is only used by step 02 (intergenic intervals, out of scope).

## 4. Module: `FILTER_TRANSCRIPT` (local)

Wraps `bin/filter_transcript.sh` (two-pass awk, preserves GFF3 gene→transcript→exon hierarchy).

```groovy
process FILTER_TRANSCRIPT {
    tag "${meta.id}.${feature_type}"
    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(gff3), val(feature_type)

    output:
    tuple val(meta), val(feature_type), path("*.${feature_type}.gff3"), emit: gff3
    path "versions.yml",                                                 emit: versions

    script:
    """
    filter_transcript.sh ${gff3} ${feature_type} > ${meta.id}.${feature_type}.gff3
    cat <<-EOF > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1)
    EOF
    """
}
```

`bin/` is auto-added to `$PATH` by Nextflow.

## 5. Module: `GFFREAD` (nf-core)

Install: `nf-core modules install gffread`

Configured in `conf/modules.config`:

```groovy
withName: 'FOMO:PREPROCESSING:GFFREAD' {
    ext.args   = '-w'
    ext.prefix = { "${meta.id}.${meta.feature_type}.exons" }
}
```

## 6. Subworkflow: `PREPROCESSING`

```groovy
workflow PREPROCESSING {
    take:
    ch_sources   // [ meta, fasta, gff3 ]

    main:
    ch_feature = ch_sources
        .combine( Channel.of('lnc_RNA', 'mRNA') )
        .map { meta, fasta, gff3, ftype ->
            tuple(meta + [feature_type: ftype], gff3, fasta, ftype)
        }

    FILTER_TRANSCRIPT(
        ch_feature.map { meta, gff3, fasta, ftype -> tuple(meta, gff3, ftype) }
    )

    ch_for_gffread = FILTER_TRANSCRIPT.out.gff3
        .join(ch_feature.map { meta, gff3, fasta, ftype -> tuple(meta, fasta) })
        .map { meta, ftype, gff3, fasta -> tuple(meta, gff3, fasta) }

    GFFREAD(
        ch_for_gffread.map { meta, gff3, fasta -> tuple(meta, gff3) },
        ch_for_gffread.map { meta, gff3, fasta -> fasta }
    )

    emit:
    spliced_fasta = GFFREAD.out.gffread_fasta
    filtered_gff3 = FILTER_TRANSCRIPT.out.gff3
    versions      = FILTER_TRANSCRIPT.out.versions.mix(GFFREAD.out.versions)
}
```

> Note: exact channel topology to be verified against the nf-core GFFREAD module signature at implementation time.

## 7. `conf/test.config`

```groovy
params {
    input  = "${projectDir}/assets/samplesheet.csv"
    outdir = "results"
}
process {
    resourceLimits = [ cpus: 2, memory: '4.GB', time: '30.min' ]
}
```

## 8. Verification

1. Run `bash bin/subsample_test_data.sh` — commit resulting `assets/test_data/`.
2. `nextflow run . -profile test,docker` (or `conda`).
3. Confirm 8 output FASTAs under `results/`: 4 sources × 2 feature types.
4. Sanity-check: `grep -c '^>' *.lnc_RNA.exons.fa` should match transcript count in matching `.lnc_RNA.gff3`.
5. Compare headers vs. legacy: run `legacy_scripts/preprocessing/00_get_feature_annotation.sh` on the same subsampled inputs. Differences expected only from the skipped longest-isoform step.

## 9. Deliberate omissions

- **Longest-isoform filtering** — skipped. Multiple isoforms per gene will appear in projection; revisit when projection results are evaluated.
- **Steps 02–05** — out of scope. Legacy scripts remain the reference.
- **AGAT statistics** — not in this PR; lands with step 05.
