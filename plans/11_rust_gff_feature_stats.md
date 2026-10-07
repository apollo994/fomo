# Plan: replace AGAT statistics with the Rust `gff-feature-stats`

> Supersedes `plans/10_replace_agat_stats.md` (the Python `bin/gff_stats_standalone.py` route).

## Context

`agat_sp_statistics.pl` is the pipeline's slowest and hungriest stat producer: **8 process
instances** across four subworkflows (`AGAT_SPSTATISTICS`, `AGAT_TARGET`, `AGAT_PROJECTED`,
`AGAT_TOP` + their `*_TO_MQC` adapters), with per-call overrides in `conf/modules.config`
reaching **12 GB / 2 h** for raw whole-genome GFFs. It also forces a gunzip of every input GFF
(`GUNZIP_RAW_SOURCE_GFF`) and a genome-FASTA join at four call sites purely to feed `--gs`.

`~/software/gff-feature-stats` (Rust, v0.2.0, `apollo994/gff-feature-stats`, single dependency
`flate2` → pure-Rust `miniz_oxide`) computes the same gene-category and transcript-type
statistics **plus `intron_stats`**, reads `.gff3`/`.gff3.gz` transparently (magic bytes, no
tabix), and runs in **~1.2 s at a flat ~20 MB RSS on a 660 MB GFF**. That removes AGAT, the
gunzip step, the FASTA joins, and every resource override in one go.

Compared with the superseded Python plan, three things change: intron columns become **real
data** instead of `NA` placeholders; the memory-batching risk disappears; and distribution
becomes a build concern instead of a `bin/`-script drop-in.

**Outcome:** richer report coverage than AGAT ever gave, no AGAT anywhere in the pipeline.

## Decisions

| Question | Decision |
|---|---|
| Distribution | **Vendor a statically linked `bin/gff-feature-stats`.** Nextflow already binds `bin/` into task containers and prepends it to `PATH` — that is how the 11 existing `bin/` scripts run today. No registry, no build step at run time. |
| `genome_coverage` | **Drop the column.** It is the only AGAT metric the tool cannot produce (needs a genome size), nothing in the pipeline reads it, and dropping it removes the last reason for FASTA plumbing near the stats steps. |
| gffread `geneID=` gap | **Leave the tool as is; document it.** `src/lib.rs:388` resolves a transcript's gene via `Parent`/`gene`/`Gene`, so `combined` and `top3` rows (gffread GTF→GFF3, no `gene` lines) report `n_genes = 0`. Transcript, exon, intron and CDS metrics are unaffected. |
| Table shape | **Redesign around the JSON**: one row per *sample × transcript type*, plus a *sample × gene category* table for input annotations. |

Explicitly **not** doing `gffread --keep-genes`: it only preserves gene records already present
in the input, and both conversions consume `gffcompare` `*.combined.gtf`, which has none — so
it would likely be a no-op while mutating two **published** deliverables and the query GFF3s
handed to `GFFCOMPARE`/`GFFCOMPARE_TOP`.

## What replaces what

| Removed | Replacement |
|---|---|
| `modules/local/agat_spstatistics.nf` | `modules/local/gff_stats.nf` → `GFF_STATS` |
| `modules/local/agat_to_mqc.nf`, `bin/agat_to_mqc.py` | `modules/local/gff_stats_to_mqc.nf` → `GFF_STATS_TO_MQC`, `bin/gff_stats_to_mqc.py` |
| aliases `AGAT_TARGET` / `AGAT_PROJECTED` / `AGAT_TOP` (+ `*_TO_MQC`) | `GFF_STATS_TARGET` / `GFF_STATS_PROJECTED` / `GFF_STATS_TOP` (+ `*_TO_MQC`) |
| `GUNZIP_RAW_SOURCE_GFF` in `subworkflows/local/preprocessing.nf` (only consumer was the AGAT channel) | nothing — the tool reads `.gz` directly |
| the four `--gs` FASTA joins + the now-unused `ch_target_fasta` params of `BENCHMARKING` and `CONSENSUS_TOP` | nothing — no genome FASTA needed |
| MultiQC sections `agat_input`, `agat_projection` | `gffstats_input`, `gffstats_genes`, `gffstats_projection` |
| the AGAT resource overrides (12 GB / 2 h, and the `gff.size() > 40_000_000` rule) | nothing — `process_single` defaults are already 300× the tool's footprint |

`bin/gff_stats_standalone.py` (currently untracked) is superseded — leave it out of git.

## 1. Build and vendor the binary

No `rustup` and no `docker` on this host (EasyBuild Rust 1.96.0-nightly, glibc 2.34,
singularity/apptainer available), so use the static-glibc route — `flate2` resolves to
pure-Rust `miniz_oxide`, so no C toolchain or musl sysroot is involved:

```sh
cd ~/software/gff-feature-stats
RUSTFLAGS='-C target-feature=+crt-static' cargo build --release
file target/release/gff-feature-stats     # expect: statically linked / static-pie
ldd target/release/gff-feature-stats      # expect: "not a dynamic executable"
cp target/release/gff-feature-stats /users/rg/fzanarello/pipelines/fomo/bin/
chmod 755 bin/gff-feature-stats
```

Cross-check it runs in the image the module will use — this is the real test of the vendoring
approach and must pass before any wiring work:

```sh
singularity exec docker://ubuntu:22.04 ./bin/gff-feature-stats -V   # → gff-feature-stats 0.2.0
```

If `crt-static` misbehaves, the fallback is a musl build inside a container
(`singularity exec docker://rust:1-alpine cargo build --release`, alpine's host target is
already musl). A dynamically linked build would *probably* work in `ubuntu:22.04` (glibc 2.34
binary on 2.35 is forward-compatible) but breaks on any older base image — don't ship that.

Provenance matters for a committed binary: record version `0.2.0` and source commit `6b31e3a`
in a comment in `modules/local/gff_stats.nf` and in the `CLAUDE.md` **Key Tools** entry, with
the rebuild recipe above. Two small upstream tidy-ups worth doing while you're there: tag
`v0.2.0` in the Rust repo, and fix `Cargo.toml`'s `repository` field (says
`fzanarello/gff-feature-stats`, the remote is `apollo994/gff-feature-stats`).

## 2. `modules/local/gff_stats.nf`

Use `ubuntu:22.04` — the container the repo's other shell-only local modules already use
(`filter_transcript.nf`, `gff_to_gene_bed.nf`, `rename_fasta_headers.nf`). The static binary
needs nothing else from the image.

```groovy
process GFF_STATS {
    tag "${meta.id}.${meta.kind}"
    label 'process_single'

    // Static build of gff-feature-stats v0.2.0 (apollo994/gff-feature-stats @ 6b31e3a),
    // vendored at bin/gff-feature-stats — Nextflow puts bin/ on PATH inside the container.
    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(gff)          // .gff3 or .gff3.gz; no genome FASTA needed

    output:
    tuple val(meta), path("*.stats.json"), emit: json
    tuple val("${task.process}"), val('gff-feature-stats'),
          eval("gff-feature-stats -V | sed 's/gff-feature-stats //'"),
          topic: versions, emit: versions_gff_feature_stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.kind}"
    """
    gff-feature-stats ${gff} ${prefix}.stats.json
    """

    stub:
    """
    echo '{"gene_category_stats":{},"transcript_type_stats":{}}' > ${task.ext.prefix ?: "${meta.id}.${meta.kind}"}.stats.json
    """
}
```

No empty-input guard needed (the AGAT module had one): a feature-less GFF yields
`{"gene_category_stats":{},"transcript_type_stats":{}}` and exit 0. The tool writes
`Saved: <path>` and any `warning:` lines to stderr — those land in `.command.err` and are
diagnostic only (see Verification §4).

## 3. `bin/gff_stats_to_mqc.py` + `modules/local/gff_stats_to_mqc.nf`

Adapter is Python (stdlib `json` only) in `container 'python:3.11'`, matching
`modules/local/filter_gff_by_id.nf`. Follow the `CLAUDE.md` routing rules: **unique filename
suffix per section, no embedded `# id:` header** (a file discovered by both `sp` and a header
renders the section twice).

CLI: `--input <json> --name <base> --role <source|target|projection> --section-id <gffstats_input|gffstats_projection> --output <tsv> [--genes-output <tsv>] [--kind <raw|filtered|decoy|projected>]`

Suffixes: `*_gffstats_input_mqc.tsv`, `*_gffstats_projection_mqc.tsv`, `*_gffstats_genes_mqc.tsv`.

**Transcript table** — one row per key in `transcript_type_stats`; `Sample` =
`<base>.<transcript_type>`, where `<base>` is the existing `ext.sample_name` value (keeps the
`1_raw`…`5_decoy_mRNA` ordering and stays unique):

| Column | Source in JSON |
|---|---|
| `role`, `kind`, `transcript_type` | CLI args / dict key |
| `n_transcripts` | `total_count` |
| `n_genes` | `associated_genes.total_count` (0 for gffread-derived files — see Decisions) |
| `n_exons`, `mean_exon_length`, `max_exon_length` | `exon_stats.total_count`, `exon_stats.length` |
| `n_introns`, `mean/min/max_intron_length` | `intron_stats.total_count`, `intron_stats.length` — **exact, not derived**; block omitted for types with no multi-exon transcript → `NA` |
| `mean/min/max_transcript_length` | `length_stats` (genomic span) |
| `mean/min/max_spliced_length` | `exon_stats.concatenated_length` (omitted for single-exon-only types → `NA`) |
| `n_cds`, `mean_cds_length` | `cds_stats` (omitted for non-coding types → `NA`) |
| `top_biotype` | highest-count key in `biotype_counts` |

`intron_stats.concatenated_length` (per-transcript intron totals) is available too — leave it
out of the table to keep it scannable, but note it in the section description.

**Gene table** (written only when `--genes-output` is given → input-side calls only) — one row
per category in `gene_category_stats`; `Sample` = `<base>.<category>`; columns `role`, `kind`,
`gene_category`, `n_genes` (`total_count`), `mean/min/max_gene_length` (`length_stats`),
`top_biotype`, `n_transcript_types` (size of `transcript_type_counts`).

Both writers emit **at least one row** (all-`NA` row, `transcript_type=none`) for an empty
JSON so MultiQC never sees a header-only file. Treat every block as optional. Log
resolved/total metric counts to stderr, as `bin/agat_to_mqc.py` does today.

Module `GFF_STATS_TO_MQC`: `label 'process_single'`, input `tuple val(meta), path(json)`,
output `tuple val(meta), path("*_mqc.tsv"), emit: tsv` (one glob catches both files), plus the
standard `topic: versions` line.

## 4. Subworkflow rewiring

**`subworkflows/local/preprocessing.nf`** — drop the `GUNZIP_RAW_SOURCE_GFF` include and call,
feed the raw source GFF3 in still-gzipped, and drop the FASTA join:

```groovy
ch_stats_gff = ch_sources.map { meta, _fa, gff3 -> tuple(meta + [kind: 'raw'], gff3) }
    .mix(ch_filtered_gff3.map { m, g -> tuple(m + [kind: 'filtered'], g) })
    .mix(RELOCATE_LOCI.out.gff3.map { m, g -> tuple(m + [kind: 'decoy'], g) })

GFF_STATS(ch_stats_gff)
GFF_STATS_TO_MQC(GFF_STATS.out.json)
```

**`subworkflows/local/projection.nf`** — replace `ch_projected_agat_in` (the `.combine()` on
`GUNZIP_TARGET.out.gunzip`) with
`ch_projected.map { meta, gff -> tuple(meta + [kind: 'projected'], gff) }`.
Keep `GUNZIP_TARGET.out.gunzip` — minimap2 and `EXTRACT_PROJECTED_LNC` still need it.

**`subworkflows/local/benchmarking.nf`** — same simplification, then **remove the
`ch_target_fasta` take** (verified: used only in the AGAT join, lines 55–61).

**`subworkflows/local/consensus_top.nf`** — same, then **remove the `ch_target_fasta` take**
(verified: used only at line 65).

**`workflows/fomo.nf`** — update both call sites:
`BENCHMARKING(ch_input.target, PROJECTION.out.gff3)` and
`CONSENSUS_TOP(ch_projected_real, BENCHMARKING.out.stats, BENCHMARKING.out.target_refs)`.
Drop `PROJECTION.out.target_fasta` from the emit list only if no consumer remains.

Update the `mqc_files` counts in each `emit:` comment — the input table now has one row per
transcript type, so the old per-file counts no longer describe it.

## 5. `conf/modules.config`

Rewrite the 8 AGAT blocks as 8 `GFF_STATS*` blocks, **preserving every `ext.prefix` and
`ext.sample_name` closure verbatim** — the `1_raw` / `2_lncRNA` / `3_decoy_lncRNA` / `4_mRNA` /
`5_decoy_mRNA` convention and the `<target>.from_<source>.<n>_<type>` projection convention
drive MultiQC's row ordering and `extra_fn_clean_exts`. Other changes:

- `ext.section_id` → `gffstats_input` / `gffstats_projection`; `ext.role` unchanged.
- Add `ext.genes_output = true` (→ `--genes-output`) on the two input-side `*_TO_MQC` blocks.
- **Delete** all AGAT `memory`/`time` overrides and the `gff.size() > 40_000_000` rule. Nothing
  replaces them: `process_single` (6 GB / 10 min) is already absurdly generous for a 20 MB-RSS,
  ~1 s tool. No `kind == 'raw'` special case.
- Delete the `FOMO:PREPROCESSING:GUNZIP_RAW_SOURCE_GFF` block.

Optionally publish the JSONs via the selector in `nextflow.config`; that alternation is
`$`-anchored, so it must read `GFF_STATS(_TARGET|_PROJECTED|_TOP)?` to catch the aliases
without also matching `GFF_STATS_TO_MQC`.

## 6. MultiQC config

`assets/multiqc/sections.yml`:
- Replace the `agat_input` / `agat_projection` `custom_data` blocks with `gffstats_input`,
  `gffstats_projection`, `gffstats_genes` (`file_format: tsv`, `plot_type: table`,
  `pconfig.id`/`namespace`). Descriptions should name `gff-feature-stats` and state the two
  caveats: introns are derived from exon gaps (`next.start − prev.end − 1`), and gffread-derived
  `combined`/`top3` rows report `n_genes = 0`.
- Replace the two `sp` patterns with the three new `fn` globs from §3.
- Update the header comment block listing the sections.

`assets/multiqc/main.yml`:
- `report_section_order`: `gffstats_input: 1000`, `gffstats_genes: 950`,
  `gffstats_projection: 900` (leave `td2_coding: 1050`, samtools and seqkit slots alone).
- `report_comment`: drop "AGAT".
- `extra_fn_clean_exts` needs no new entry (`_mqc`, `.tsv` are there) — but confirm the new
  `.<transcript_type>` sample-name suffixes are not accidentally stripped.

## 7. Deletions & docs

Delete `modules/local/agat_spstatistics.nf`, `modules/local/agat_to_mqc.nf`,
`bin/agat_to_mqc.py`. `grep -rniI agat --include='*.nf' --include='*.config' --include='*.yml'
--include='*.json' .` should then hit only `legacy_scripts/` (reference implementations — leave
them alone).

`CLAUDE.md`: **Key Tools** (AGAT → `gff-feature-stats`, with the vendored-binary provenance and
rebuild recipe; keep AGAT named only as the legacy reference), the **Subworkflow map** table
("AGAT + SeqKit stats" → "GFF stats + SeqKit"), the **Preprocessing detail** statistics row, and
the **Reporting conventions** paragraph (the MultiQC-container-reuse note no longer applies to
these adapters; `modules/nf-core/multiqc/` and `modules/local/samtools_to_mqc.nf` remain the
version-bump sites).

## Verification

1. **Upstream tests**: `cd ~/software/gff-feature-stats && cargo fmt --all -- --check &&
   cargo clippy --all-targets -- -D warnings && cargo test --all` — the integration test
   compares against the pinned `tests/data/sample.expected.json`.
2. **Static-binary pre-flight** (§1): `file`/`ldd` show static, and
   `singularity exec docker://ubuntu:22.04 ./bin/gff-feature-stats -V` prints the version.
   Everything downstream depends on this.
3. **Numeric parity vs AGAT, before deleting anything**: from a previous full run's
   `work/`/`results/`, take 3–4 GFFs spanning raw source, filtered lncRNA, decoy and a projected
   model; run the binary and diff against the matching `*.stats.yaml`.
   - Must match: transcript counts, exon counts, transcript/exon length extrema.
   - **Now checkable (new):** introns against AGAT's *number of introns in exon*, *shortest /
     longest intron*, and *total intron length per exon* — this is the metric that had no
     cross-check in the Python plan.
   - Expected deltas to record, not fix: `n_genes` can be lower than AGAT's (genes whose
     transcripts have no exon/CDS records never enter `gene_category_stats`), and gffread-derived
     files report `n_genes = 0`.
4. **Warning audit**: `gff-feature-stats` emits a `warning:` line per transcript whose declared
   span differs from its exon envelope. `BAM_TO_GFF` sets each transcript's span to exactly the
   first-to-last exon block, so projected models should produce **zero** warnings — any warning
   there is a genuine `bin/bam_to_gff.sh` signal worth chasing. Warnings on raw Ensembl input are
   informational.
5. **Wiring**: `nextflow run . -profile test,docker -stub` — completes with no
   `AGAT`/`GUNZIP_RAW_SOURCE_GFF` tasks and no channel-arity errors from the removed
   `ch_target_fasta` params.
6. **Full smoke run**: `nextflow run . -profile test,docker`, then inspect
   `results/multiqc/multiqc_report.html`:
   - the three new sections each render **exactly once**;
   - input table has one row per (species × kind × transcript type) with the `1_raw`…
     `5_decoy_mRNA` ordering intact;
   - projection table covers the 20 per-source + combined models and the 2 top3 rows;
   - intron columns are **populated** (not `NA`) for every multi-exon type — this is the headline
     difference from the previous plan;
   - `results/multiqc/multiqc_data/` contains the parsed `gffstats_*` tables.
7. **Confirm the motivation**: compare `realtime`/`peak_rss` for `GFF_STATS` in
   `results/pipeline_info/execution_trace.txt` against the AGAT rows from an earlier run.
8. `grep -c AGAT` over `.nf`/`.config`/`.yml` outside `legacy_scripts/` → 0.

## Risks / follow-ups

- **Vendored binary is linux-x86_64 only.** Every container-based profile is fine; a bare-metal
  macOS/arm run would need a locally rebuilt binary. Document the recipe (§1) next to the
  `CLAUDE.md` Key Tools entry. Revisit if the tool ever lands on bioconda/crates.io — then this
  becomes a one-line `container` pin and the binary leaves git.
- **`n_genes = 0` on `combined`/`top3` rows** by decision. If it starts to matter, the fix is
  ~2 lines at `src/lib.rs:388` (add `geneID`/`gene_id` to the fallback chain) released as v0.3.0
  — `intron_stats` already made the output a deliberate superset of annotrieve, so there is no
  byte-compat constraint left to protect.
- **Projected genes all land in `non_coding`**: `BAM_TO_GFF` emits `gene`/`transcript`/`exon`
  with no CDS. Expected — worth one line in the section description so the gene table is not
  misread.
- **Intron semantics**: gaps are used as-is, so overlapping exons yield short or negative gaps
  rather than an error. Pipeline inputs are Ensembl-style or generated by our own tools, so this
  is a theoretical concern; §4's warning audit is the tripwire.
- `genome_coverage` is gone by decision. Restoring it would mean a genome-size argument plus a
  merged-interval pass in Rust, and reinstating the FASTA join at all four call sites.
