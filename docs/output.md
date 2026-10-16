# Outputs

Everything is published under `--outdir`. `<track>` is a feature type with an optional decoy
suffix: `lnc_RNA`, `lnc_RNA.decoy`, `mRNA`, `mRNA.decoy` (mRNA with `--include_mrna`, decoys with
`--include_decoy`).

```
<outdir>/
├── targets/<target>/
│   ├── alignment/   <target>.allModels.<track>.raw.bam  (+ .bam.csi)
│   ├── annotation/  <target>.allModels.<track>.raw.gff3
│   │                <target>.from_<source>.<track>.projected.gff3
│   ├── gffcompare/  <target>.from_<model>.<track>.gffcompare.{stats,tracking,loci,annotated.gtf}
│   │                <target>.from_<model>.<track>.gffcompare.<query>.{tmap,refmap}
│   ├── curated/     <target>.curated.<track>.gff3.gz
│   │                <target>.curated.<track>.merged.gff3.gz
│   │                <target>.curated.<track>.curation.{tsv,json}
│   └── multiqc/     multiqc_report.html, multiqc_report_data/
├── sources/<source>/
│   └── annotation/  <source>.<feature type>.filtered.gff3
├── summary/
│   ├── multiqc/     fomo_run_summary.html, fomo_run_summary_data/
│   └── tables/      run_*_mqc.tsv, run_*_mqc.json, run_summary.json
└── pipeline_info/   execution_{report,timeline,trace}_<time>.*, pipeline_dag_<time>.html
```

## targets/&lt;target&gt;/

**alignment/** — one coordinate-sorted BAM per track: every source's transcripts aligned to the
target in one minimap2 run (`-x splice -uf --cs`). A CSI index, not BAI, so chromosomes longer
than 2²⁹ bp work.

**annotation/** — the projected models.
- `allModels.<track>.raw.gff3` pools every source except the target itself. For lncRNA the models
  TD2 calls coding are already removed.
- `from_<source>.<track>.projected.gff3` holds one source's share of the same models.

Rows are `gene` → `transcript` → `exon`. A transcript's ID is `<source transcript>|<source>`,
unique across sources, and carries `source`, `gene_class`, `query_length`, `aligned_bases`,
`pct_aligned`, `projected_length`, `pct_projected` and minimap2's `NM`, `AS` and `de`. Single-exon
alignments are left out unless `--include_single_exon`.

**gffcompare/** — for targets of role `both` only: each source's models (`<model>` = the source)
and the pool (`<model>` = `allModels_raw`) against the target's filtered annotation
(`gffcompare -M --no-merge`). `.stats` holds sensitivity and precision per level; decoy tracks are
compared with the real reference, so their precision is the false-positive baseline.

**curated/** — the curated annotation: intron chains shared exactly by at least
`--curate_min_species` source species, minus models with an exon on a same-strand reference exon.
- `curated.<track>.gff3.gz` has one gene per group of chains that share a splice junction and one
  transcript per chain.
  - Gene attributes: `ID`, `biotype`, `n_chains`, `n_species`, `ref_location` (`antisense_exonic`,
    `intronic`, `sense_span_overlap`, `antisense_intronic`, `divergent`, `intergenic`,
    `no_reference`) and the reference gene it refers to (`ref_gene_id`, `ref_gene_name`,
    `ref_gene_biotype`, `ref_gene_orientation`, `ref_distance`).
  - Transcript attributes: `representative_transcript`, `supporting_transcripts`,
    `supporting_species`, `n_sj`.
  - IDs are `fomo_<tag>_g_<hash>`, `fomo_<tag>_t<n>_<hash>`, `fomo_<tag>_t<n>e<k>_<hash>`, where the
    hash is computed from the gene's intron chains: the same chains on the same assembly always
    get the same IDs, whatever the source set.
- `curated.<track>.merged.gff3.gz` is the target's reference with the curated genes inserted
  after each seqid's last reference line. It exists only for targets with a gff3.
- `curated.<track>.curation.{tsv,json}` report what was observed, what the reference filter
  removed and what was added.

**multiqc/** — the target's report: annotation statistics of every source (filtering funnel, TD2,
GFF statistics, FASTA statistics), the projection (model counts, alignment), the accuracy
scatter, the curation and the software versions.

## sources/&lt;source&gt;/annotation/

`<source>.<feature type>.filtered.gff3` lists exactly the transcripts that were projected:
stranded, spliced, the longest isoform per gene (unless `--longest_isoform false`) and, for lncRNA,
not coding according to TD2. For a `both` species it is also the benchmarking reference.

## summary/

The run-level report, aggregated across targets and sources: species per role, the source
filtering funnel, genes per target before and after curation, the source × target accuracy
heatmaps, pooled models against the best single donor, and a donor ranking. `tables/` holds the
same numbers as TSV/JSON; `run_summary.json` holds all of them in one file.

## pipeline_info/

Nextflow's report, timeline, trace and DAG. Each run gets its own files, named with the launch
time (`--trace_report_suffix` sets another suffix).
