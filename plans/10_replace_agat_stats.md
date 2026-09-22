# Plan: replace AGAT statistics with `bin/gff_stats_standalone.py`

## Context

`agat_sp_statistics.pl` is the pipeline's slowest and hungriest stat producer. It runs
**8 process instances** across four subworkflows (`AGAT_SPSTATISTICS`, `AGAT_TARGET`,
`AGAT_PROJECTED`, `AGAT_TOP` + their `*_TO_MQC` adapters) and `conf/modules.config`
carries per-call resource overrides up to **12 GB / 2 h** for raw whole-genome GFFs.
It also forces a gunzip of every input GFF (`GUNZIP_RAW_SOURCE_GFF`) and a genome-FASTA
join at four call sites purely to feed `--gs`.

`bin/gff_stats_standalone.py` computes equivalent per-transcript-type and per-gene-category
statistics in one dependency-free sequential pass, reads `.gff3` and `.gff3.gz` transparently
(gzip detected by magic bytes), and emits JSON. Replacing AGAT removes the AGAT container,
the gunzip step, the FASTA joins, and the heavy resource overrides.

**Outcome:** identical report coverage (richer, in fact), no AGAT anywhere in the pipeline.

Prerequisite: `git add bin/gff_stats_standalone.py` — currently untracked (already mode 755).

## Decisions

| Question | Decision |
|---|---|
| Metrics the JSON cannot supply (intron lengths, `genome_coverage`) | **Keep the columns, render `NA`.** No new computation added to `gff_stats_standalone.py`. `n_introns` is derived in the adapter (`n_exons − n_transcripts`) and flagged as approximate. |
| Table shape | **Redesign around the JSON.** One row per *sample × transcript type* (transcript tables) plus a *sample × gene category* table for input annotations. |
| gffread consensus GFF3s (no `gene` lines, gene in `geneID=`) | **Script-side attribute fallback, gated behind a flag** — *not* `gffread --keep-genes`. Rationale below. |

### Why not `gffread --keep-genes`

`--keep-genes` only preserves gene records **present in the input**. Both conversions
(`COMBINED_GTF_TO_GFF`, `COMBINED_TOP_GTF_TO_GFF`) consume `gffcompare` `*.combined.gtf`,
which contains only `transcript`/`exon` records — there are no gene records to keep, so the
flag would likely be a no-op. It would also alter two **published** deliverables
(`nextflow.config` publishes both conversions) and change the query GFF3s handed to
`GFFCOMPARE`/`GFFCOMPARE_TOP`. The script-side fallback has zero downstream effect: it only
changes what the stats reader sees.

## What replaces what

| Removed | Replacement |
|---|---|
| `modules/local/agat_spstatistics.nf` | `modules/local/gff_stats.nf` → `GFF_STATS` |
| `modules/local/agat_to_mqc.nf`, `bin/agat_to_mqc.py` | `modules/local/gff_stats_to_mqc.nf` → `GFF_STATS_TO_MQC`, `bin/gff_stats_to_mqc.py` |
| aliases `AGAT_TARGET` / `AGAT_PROJECTED` / `AGAT_TOP` (+ `*_TO_MQC`) | `GFF_STATS_TARGET` / `GFF_STATS_PROJECTED` / `GFF_STATS_TOP` (+ `*_TO_MQC`) |
| `GUNZIP_RAW_SOURCE_GFF` in `subworkflows/local/preprocessing.nf` (its only consumer was the AGAT channel) | nothing — the script reads `.gz` directly |
| the four `--gs` FASTA joins, and the now-unused `ch_target_fasta` params of `BENCHMARKING` and `CONSENSUS_TOP` | nothing — no genome FASTA needed |
| MultiQC sections `agat_input`, `agat_projection` | `gffstats_input`, `gffstats_genes`, `gffstats_projection` |

## 1. `bin/gff_stats_standalone.py` — one gated addition

Add an opt-in CLI flag (e.g. `--gene-attr-fallback`) that extends the transcript→gene
lookup in `compute_features_statistics` (currently line 256):

```python
# default (unchanged, byte-identical to annotrieve):
attr.get('Parent') or attr.get('gene') or attr.get('Gene') or None
# with --gene-attr-fallback:
... or attr.get('geneID') or attr.get('gene_id') or attr.get('locus') or None
```

Default off keeps the annotrieve byte-for-byte guarantee intact; `GFF_STATS` always passes
the flag. Effect: consensus/top3 GFF3s report real associated-gene counts instead of 0.
Note the residual limitation in a comment — with no `gene` feature line those gene IDs never
enter `gene_info`, so `gene_category_stats` stays empty for gffread-derived files and the
transcript table's `n_genes` comes from `associated_genes.total_count`.

## 2. `modules/local/gff_stats.nf`

Mirror `modules/local/relocate_loci.nf` (the existing `python:3.11` local-module pattern):

```groovy
process GFF_STATS {
    tag "${meta.id}.${meta.kind}"
    label 'process_single'
    container 'python:3.11'

    input:
    tuple val(meta), path(gff)          // NOTE: no genome FASTA, .gz accepted

    output:
    tuple val(meta), path("*.stats.json"), emit: json
    tuple val("${task.process}"), val('python'),
          eval('python3 --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '--gene-attr-fallback'
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    gff_stats_standalone.py ${gff} ${prefix}.stats.json ${args}
    """

    stub:
    """
    echo '{"gene_category_stats":{},"transcript_type_stats":{}}' > ${task.ext.prefix ?: "${meta.id}.${meta.kind}"}.stats.json
    """
}
```

The AGAT module's empty-GFF guard is unnecessary — the script returns
`{"gene_category_stats":{},"transcript_type_stats":{}}` on a feature-less file. Switch the
script's `main()` from raw `sys.argv` indexing to `argparse` so the flag can be passed
alongside the two positional paths in any order.

## 3. `bin/gff_stats_to_mqc.py` + `modules/local/gff_stats_to_mqc.nf`

Adapter reads the JSON (stdlib `json` only → `python:3.11` container, no MultiQC image
needed) and writes MultiQC custom-content TSVs. Follow the routing rules in `CLAUDE.md`:
**unique filename suffix per section, no embedded `# id:` header.**

CLI: `--input <json> --name <base> --role <source|target|projection> --section-id <gffstats_input|gffstats_projection> --output <tsv> [--genes-output <tsv>] [--kind <raw|filtered|decoy|projected>]`

Suffixes: `*_gffstats_input_mqc.tsv`, `*_gffstats_projection_mqc.tsv`, `*_gffstats_genes_mqc.tsv`.

**Transcript table** — one row per transcript type present in `transcript_type_stats`;
row `Sample` = `<base>.<transcript_type>` (keeps the existing numeric ordering prefix from
`ext.sample_name` and stays unique):

| Column | Source in JSON |
|---|---|
| `role`, `kind`, `transcript_type` | CLI args / dict key |
| `n_transcripts` | `total_count` |
| `n_genes` | `associated_genes.total_count` |
| `n_exons` | `exon_stats.total_count` |
| `n_introns` | derived `max(0, n_exons − n_transcripts)` — approximate; a CDS-only transcript with no `exon` records makes it an undercount |
| `mean/min/max_transcript_length` | `length_stats` (genomic span) |
| `mean/max_exon_length` | `exon_stats.length` |
| `mean/min/max_spliced_length` | `exon_stats.concatenated_length` (absent for single-exon-only types → `NA`) |
| `n_cds`, `mean_cds_length` | `cds_stats` (absent for non-coding types → `NA`) |
| `mean/min/max_intron_length`, `genome_coverage` | **always `NA`** — placeholders, see Decisions |
| `top_biotype` | highest-count key in `biotype_counts` |

**Gene table** (written only when `--genes-output` is given → input-side calls only) — one row
per category in `gene_category_stats`; `Sample` = `<base>.<category>`; columns `role`, `kind`,
`gene_category`, `n_genes` (`total_count`), `mean/min/max_gene_length` (`length_stats`),
`top_biotype`, `n_transcript_types` (len of `transcript_type_counts`).

Both writers must emit **at least one row** (all-`NA` placeholder row with
`transcript_type=none`) for an empty JSON, so MultiQC never sees a header-only file.
Log resolved/total metric counts to stderr like `bin/agat_to_mqc.py` does today.

Module `GFF_STATS_TO_MQC`: `label 'process_single'`, `container 'python:3.11'`,
input `tuple val(meta), path(json)`, output `tuple val(meta), path("*_mqc.tsv"), emit: tsv`
(glob catches both files), plus the standard `topic: versions` line.

## 4. Subworkflow rewiring

**`subworkflows/local/preprocessing.nf`** — drop the `GUNZIP_RAW_SOURCE_GFF` include and
call; feed the raw source GFF3 (still `.gz`) straight in, and drop the FASTA join:

```groovy
ch_stats_gff = ch_sources.map { meta, _fa, gff3 -> tuple(meta + [kind: 'raw'], gff3) }
    .mix(ch_filtered_gff3.map { m, g -> tuple(m + [kind: 'filtered'], g) })
    .mix(RELOCATE_LOCI.out.gff3.map { m, g -> tuple(m + [kind: 'decoy'], g) })

GFF_STATS(ch_stats_gff)
GFF_STATS_TO_MQC(GFF_STATS.out.json)
```

**`subworkflows/local/projection.nf`** — replace `ch_projected_agat_in` (the
`.combine()` on `GUNZIP_TARGET.out.gunzip`) with
`ch_projected.map { meta, gff -> tuple(meta + [kind: 'projected'], gff) }`.
`GUNZIP_TARGET.out.gunzip` is still needed for minimap2 / `EXTRACT_PROJECTED_LNC`, so keep it.

**`subworkflows/local/benchmarking.nf`** — same simplification; then **remove the
`ch_target_fasta` take** (verified: used only in the AGAT join at lines 55–61).

**`subworkflows/local/consensus_top.nf`** — same; **remove the `ch_target_fasta` take**
(verified: used only at line 65).

**`workflows/fomo.nf`** — update both call sites:
`BENCHMARKING(ch_input.target, PROJECTION.out.gff3)` and
`CONSENSUS_TOP(ch_projected_real, BENCHMARKING.out.stats, BENCHMARKING.out.target_refs)`.
Keep `PROJECTION.out.target_fasta` in the emit list only if another consumer remains —
otherwise delete that emit too.

Update the `mqc_files` comments/counts in each `emit:` block (row counts change: the input
transcript table now has one row per type).

## 5. `conf/modules.config`

Rewrite the 8 AGAT blocks as 8 `GFF_STATS*` blocks, **preserving every existing
`ext.prefix` and `ext.sample_name` closure verbatim** (the `1_raw` / `2_lncRNA` /
`3_decoy_lncRNA` / `4_mRNA` / `5_decoy_mRNA` ordering convention and the
`<target>.from_<source>.<n>_<type>` projection convention must survive — `sections.yml`
sort order and `extra_fn_clean_exts` depend on them). Changes beyond the rename:

- `ext.section_id` → `gffstats_input` / `gffstats_projection`.
- `ext.role` unchanged (`source` / `target` / `projection`).
- Add `ext.genes_output = true` (or set `--genes-output` via `ext.args`) on the two
  input-side `*_TO_MQC` blocks only.
- **Delete** the AGAT resource overrides (`memory = { meta.kind == 'raw' ? 12.GB ... }`,
  `time`, and the `gff.size() > 40_000_000` rule on `AGAT_PROJECTED`). Start from
  `process_single` defaults (6 GB / 10 min) and add a single modest override for
  `kind == 'raw'` **only if** the benchmark in Verification §2 shows it's needed.
- Delete the `FOMO:PREPROCESSING:GUNZIP_RAW_SOURCE_GFF` block.

Optionally add the JSONs to the publish selector in `nextflow.config`; note the alternation
is `$`-anchored, so it must read `GFF_STATS(_TARGET|_PROJECTED|_TOP)?` to catch the aliases
without also matching `GFF_STATS_TO_MQC`.

## 6. MultiQC config

`assets/multiqc/sections.yml`:
- Replace the `agat_input` / `agat_projection` `custom_data` blocks with `gffstats_input`,
  `gffstats_projection`, `gffstats_genes` (`file_format: tsv`, `plot_type: table`,
  `pconfig.id`/`namespace`, descriptions naming `bin/gff_stats_standalone.py` instead of AGAT).
- Replace the two `sp` patterns with the three new `fn` globs from §3.
- Update the header comment block listing the sections.

`assets/multiqc/main.yml`:
- `report_section_order`: `gffstats_input: 1000`, `gffstats_genes: 950`,
  `gffstats_projection: 900` (keep `td2_coding: 1050`, samtools/seqkit slots unchanged).
- `report_comment`: drop "AGAT".
- `extra_fn_clean_exts` needs no new entry (`_mqc`, `.tsv` already listed) — but confirm the
  new `.<transcript_type>` sample suffixes are not accidentally stripped.

## 7. Deletions & docs

Delete `modules/local/agat_spstatistics.nf`, `modules/local/agat_to_mqc.nf`,
`bin/agat_to_mqc.py`. Grep for stragglers: `grep -rniI agat --include='*.nf'
--include='*.config' --include='*.yml' --include='*.json' .` should return only
`legacy_scripts/` (reference implementations — leave untouched).

`CLAUDE.md`: update the **Key Tools** list (AGAT → `bin/gff_stats_standalone.py`, keep AGAT
mentioned only as the legacy reference), the **Subworkflow map** table ("AGAT + SeqKit stats"
→ "GFF stats + SeqKit"), the **Preprocessing detail** table's statistics row, and the
**Reporting conventions** paragraph (the MultiQC-container reuse note no longer applies to the
new adapters — they use `python:3.11`; `modules/nf-core/multiqc/` and
`modules/local/samtools_to_mqc.nf` remain the version-bump sites).

## Verification

1. **Numeric parity check** (do this first, before deleting AGAT): pick an existing results
   dir / work dir from a previous full run and, for 3–4 GFFs spanning raw source, filtered
   lncRNA, decoy, and a projected model, run `bin/gff_stats_standalone.py` and diff against
   the matching `*.stats.yaml`. Confirm `n_transcripts`, `n_exons`, transcript/exon lengths
   match, and **document the expected deltas**: genes whose transcripts have no exon/CDS
   records are excluded from `gene_category_stats` (AGAT counts them), so `n_genes` can be
   lower; gffread-derived files report genes only via `associated_genes`.
2. **Benchmark the motivation**: time + peak RSS for the script vs the AGAT task on the
   largest raw source annotation available, then compare against `memory`/`realtime` in
   `results/pipeline_info/execution_trace.txt` from a previous run. Use the result to set (or
   justify omitting) the `kind == 'raw'` resource override.
3. **Wiring**: `nextflow run . -profile test,docker -stub` — must complete with no
   `AGAT`/`GUNZIP_RAW_SOURCE_GFF` tasks and no channel-arity errors from the removed
   `ch_target_fasta` params.
4. **Full smoke run**: `nextflow run . -profile test,docker`. Then inspect
   `results/multiqc/multiqc_report.html`:
   - three new sections present, each rendered **exactly once** (the double-render trap in
     `CLAUDE.md`);
   - input transcript table has one row per (species × kind × transcript type) with the
     `1_raw`…`5_decoy_mRNA` ordering preserved;
   - projection table covers all 20 per-source + combined models plus the 2 top3 rows;
   - no row is entirely `NA`; the intron-length and `genome_coverage` columns are uniformly
     `NA` as intended;
   - `results/multiqc/multiqc_data/` contains the parsed `gffstats_*` tables.
5. `grep -c AGAT` over `.nf`/`.config`/`.yml` outside `legacy_scripts/` → 0.

## Risks / follow-ups

- **Batching memory**: the script only flushes at a seqid boundary *after* 200 k lines, so a
  single very large first chromosome is held entirely in memory. Fine for the test data and
  for insect/vertebrate Ensembl GFFs; revisit if a pathological single-sequence annotation
  shows up.
- **Projected-model gene categories**: `BAM_TO_GFF` emits `gene`/`transcript`/`exon` with no
  CDS, so every projected gene lands in `non_coding` — expected, worth a one-line note in the
  section description so the report isn't misread.
- Intron lengths and `genome_coverage` stay `NA` by decision; if they're wanted later, add
  them as an additive opt-in block in the script (intron lengths from exon coordinates,
  coverage from merged gene intervals + a `.fai`) without touching the default output.
