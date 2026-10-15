# How FOMO works

`workflows/fomo.nf` validates the samplesheet and splits it into sources (`source`, `both`) and
targets (`target`, `both`) with two filters, so a `both` row reaches both. It then runs six
subworkflows; each emits `data` (statistics records) and, where a native MultiQC module reads
them, `mqc` (tool files).

| subworkflow | per | does |
|:---|:---|:---|
| PREPROCESSING | source | candidate transcripts, spliced FASTA, TD2, the filtered annotation, decoys |
| PROJECTION | target × track | one alignment of every source, BAM → GFF3, TD2 on the pooled lncRNA, split per source |
| BENCHMARKING | `both` target × track | gffcompare of every model set against the target's filtered annotation |
| CURATION | target × track | chains shared by several species, filtered by and merged into the reference |
| REPORTING | target | one MultiQC report |
| RUN_SUMMARY | run | one aggregate report |

Preprocessing is target-independent (its cost grows with S, not S·T); everything per target runs
after the sources are merged.

## Tracks and meta

A track is a feature type plus the decoy flag. It is stamped on the meta map once
(`withTrack`, `subworkflows/local/fomo_utils.nf`): `track` (`lnc_RNA`, `lnc_RNA.decoy`, `mRNA`,
`mRNA.decoy`, the file-name segment), `label` (`lncRNA`, `decoy_lncRNA`, ...), `order` (row
order in reports) and `td2` (whether TD2 filters it). The enabled feature types are derived once,
in `workflows/fomo.nf`.

Other meta keys: `id` (the species; for target-scoped tasks the target), `target_id` (set on
everything that belongs to one target: it routes MultiQC content and names the publish
directory), `role`, `kind` (`raw`, `filtered`, `decoy`, `projected`) on annotation statistics,
`stage` (`input`, `projected`) on TD2, and `exclude_self` on the alignment.

Channels are paired with `pairBy`/`pairById` on an explicit key. The projection ⋈ reference
pairing keys on `[target, feature type]`: on the feature type alone a projection would be
compared with every target's reference.

## Preprocessing

`bin/filter_transcript.py` keeps the stranded, spliced transcripts of a feature type and, with
`--longest_isoform`, the longest spliced isoform per gene; it is the only place transcripts are
selected. gffread extracts their sequences; `RENAME_FASTA_HEADERS` names each record
`>tid|type|species` and writes the map back to the GFF3 ID; TD2 drops lncRNA with a complete ORF
whose PSAURON score reaches `--td2_psauron_min`. `FILTER_ANNOTATION` then cuts the candidate GFF3
to exactly the records of that final FASTA, matched through the map by exact ID, and fails if the
three disagree. That filtered annotation is published, is the benchmarking reference of a `both`
species, and describes exactly what was projected.

Decoys (`--include_decoy`) move each candidate model's shape into the source's own intergenic
space: the complement of every gene-level feature of the raw annotation, coding genes included,
whatever tracks are enabled. Placement is random but seeded (`--relocate_seed`).

## Inputs stay compressed

No task writes a decompressed genome or annotation to the shared work directory. gff-feature-stats,
`filter_transcript.py` and `curate_models.py` read `.gz` directly; minimap2 indexes a gzipped
genome. gffread cannot read gzip and reads its genome twice, so `EXTRACT_TRANSCRIPTS` stages the
genome once, plain, on node-local scratch (`$TMPDIR`). `samtools faidx` is not used: it accepts
bgzip but rejects the plain gzip that Ensembl and NCBI distribute.

## Projection: merge, align once, split

Every source's spliced FASTA of a track is concatenated into one file and aligned to each target
in one minimap2 run. The species travels with each record: `>tid|type|species` becomes the BAM
read name, `bin/bam_to_gff.sh` writes it as `source=<species>` on the gene and transcript rows and
names the transcript `<tid>|<species>` (transcript IDs repeat across species, and every source
shares one GFF3), and `bin/gff_by_source.py` splits the result back per source. The split also
writes a manifest (`file`, `source`, the pooled file first) from which the batch tasks take each
file's identity; no step parses a file name for it.

A species is never projected onto itself: for a `both` target `MINIMAP2_ALIGN` removes the
target's own records from the query (`bin/fasta_drop_species.sh`, reading the same species field).
A target onto which nothing projects is reported with zero models rather than dropped.

## Statistics records and reports

Every statistics producer writes data records (`fomo-data/1`, [data_records.md](data_records.md))
in the same task. The per-target MULTIQC task renders its custom content from records
(`bin/render_mqc.py`, sections registered in `bin/fomo_sections.py`); RUN_SUMMARY_TABLES builds the
run summary from the same records. Records with a target identity go to that target's report,
source-side records to every report.

## Containers, configuration and versions

Each module carries one `img_*` label; `conf/containers.config` maps labels to pinned images.
`conf/base.config` holds every resource request, `conf/modules.config` every `ext.*` and
`publishDir` (selectors are `.*:NAME$`), and site profiles only site settings. Optional branches
(decoys, curation) are gated by giving them an empty input channel. Every module puts its tool
versions on the `versions` topic; `fomo.nf` collects them into the Software Versions section.
