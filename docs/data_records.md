# Data records

Every statistics producer writes its numbers as records, in the task that measured them. The
reports are rendered from records only. A record file is `*.data.json` and holds one record or a
list:

```json
{"schema": "fomo-data/1",
 "section": "gffstats",
 "identity": {"species": "Lycaena_alciphron_282377", "side": "target", "target": "Lycaena_alciphron_282377",
              "source": "Lycaena_thersamon_265360", "track": "lnc_RNA", "label": "lncRNA", "order": 1,
              "kind": "projected"},
 "data": {}}
```

Identity keys: `species`, `side` (`source` or `target`), `target`, `source`, `track`, `label`,
`order`, `kind`, `stage`. `bin/fomo_stats.py` defines the sections and the keys each requires
(`SECTIONS`) and validates every record it reads; `tests/test_docs.py` keeps this page in step.

### transcript_filter

FILTER_TRANSCRIPT, per source and feature type. Requires `species`, `side`, `track`.
Data: `n_transcripts`, `n_stranded`, `n_spliced`, `n_selected`.

### gffstats

GFF_STATS on every source's raw, filtered and decoy annotation, on a pure target's own annotation,
and on every projected model set (with `source`). Requires `species`, `side`, `kind`.
Data: the gff-feature-stats JSON, verbatim.

### td2_coding

TD2_CODING_FILTER, per source (`stage` input) and per target (`stage` projected). Requires
`species`, `side`, `track`, `stage`. Data: `n_in`, `n_coding`, `n_kept`, `psauron_min`.

### alignment

SAMTOOLS_STATS, per target and track. Requires `target`, `side`, `track`. Data:
`raw_total_sequences`, `reads_mapped`, `reads_mapped_percent`, `non_primary_alignments`,
`error_rate`, `average_quality`, `average_length`.

### fasta_stats

SEQKIT_STATS, per source and track. Requires `species`, `side`, `track`. Data: every column of
`seqkit stats --all` but the file name.

### accuracy

GFFCOMPARE_BATCH, per target, track and model (`source` is a species or `allModels_raw`).
Requires `target`, `source`, `side`, `track`. Data: `levels` (sensitivity `sn` and precision
`pr` per level: Base, Exon, Intron, Intron_chain, Transcript, Locus) and `counts` (query and
reference transcripts and loci, matching intron chains, transcripts and loci, missed and novel
exons, introns and loci with their totals).

### curation

CURATE_MODELS, per target and track. Requires `target`, `side`, `track`. Data: `observed`,
`discarded_by_reference`, `novel` (as in the published `curation.json`) and `has_reference`.
