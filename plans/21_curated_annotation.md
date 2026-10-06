# 21 — Curated annotation: chains shared by ≥ 2 species, merged with the reference

## Status: phase A done (2026-10-06) on branch `curated-annotation` (uncommitted); phases B–C pending

Phase A = `bin/curate_models.py` + `tests/curate_models/test_curate_models.py` (43 fixture checks
pass). The real-data validation is in **Phase A results** at the end.

## Decisions (2026-10-06)

| Question | Decision |
|---|---|
| What can be reported | **Only intron chains shared exactly (every SJ identical) by ≥ `curate_min_species` (default 2) species.** Other models are ignored for reporting. |
| SJ | One intron, `(chr, strand, donor, acceptor)`, exact coordinates, no tolerance. Single-exon models have no SJ and are dropped. |
| Grouping into genes | Supported chains that **share ≥ 1 SJ** form one gene (connected components). **Only supported chains link**: unsupported models never bridge two chains. |
| Representative (one per gene) | Over every model of every supported chain in the gene: **most SJs → the chain shared by the most species (after the reference filter) → longest exonic length → lowest `de` → ID**, so the most complete supported model is reported. The later keys also pick which model (which ends) represents the chain. *Rev 2 briefly put species first; reverted the same day (rev 3) at the user's request.* |
| Reference filter | Applied **per model, before support is counted**. A model is removed if any of its exons overlaps **any exon entry of the reference on the same strand**, whatever the parent type: mRNA, lncRNA, antisense_RNA, misc_RNA `transcript`, tRNA, snoRNA, pseudogenic, … A transcript-like feature with no exon children counts by its span. The overlapped parent types are recorded for the report. `--strand both` exists for standalone use. |
| Self projections | **Never count as support.** They are rare after the filter anyway, since they mostly overlap the reference models they came from. Same reasoning as the ranking-pool filter in `consensus_top.nf`. CLI flag `--keep-self` for standalone use only. |
| Tracks | **Every** `allModels.<gtype>.raw.gff3`: lncRNA, mRNA (`--include_mrna`), and their decoys (`--include_decoy`), i.e. F·D·T tasks. `top3` is not curated. |
| Deliverable | Per track, **gzipped**: `<target>.curated.<gtype>.gff3.gz` (curated genes only) **and** `<target>.curated.<gtype>.merged.gff3.gz` (target reference + curated genes) when the target has a reference. |
| Biotype labels | Prefixed with `fomo_` on gene (`gene_biotype`) and transcript (`transcript_biotype`): `fomo_lncRNA`, `fomo_mRNA`, `fomo_decoy_lncRNA`, `fomo_decoy_mRNA`. |
| Report | Per track, three tables (§4): **(1)** multi-species chains/genes observed (before the reference filter); **(2)** how many were discarded by overlap with the reference, broken down by the overlapped reference type; **(3)** summary stats of the novel (kept) models. Phase A writes them as TSV; the pipeline turns them into MultiQC sections. |
| Representative attributes | `n_supporting_species`, `supporting_species`, `supporting_models` (IDs of every model with the representative's exact chain; IDs carry the species, `<tid>\|<source>`), plus `n_sj`, `n_chains` / `n_species_gene` on the gene. |
| Location vs reference (rev 2) | `ref_location=<class>` on gene **and** transcript; `ref_nearest` / `ref_distance` on the transcript (§2a). |
| Tools | **gffcompare** decides exact chain identity (multi-input tracking, `--no-merge`). **Python** handles ref overlap, SJ graph, representative and output. **gffread** validates the curated and merged GFF3s. One process, one container. |

## Context

`tools/fomo-curate/curate_shared_chains.py` (outside the pipeline) curates a target's raw
models by grouping intron chains within ±10 bp. The review of 2026-10-05 ran it on 82 targets
(all acari/apiaceae/solanaceae/drosophilidae, plus 3 lepidoptera and 3 mollusca targets) against
gffcompare (`tools/fomo-curate/gffcompare_shared_chains.py`):

- **At tol 0 it is identical to gffcompare** on all 82 targets: same loci, members and
  representatives. There are no fake exact chains.
- **At tol 10, the default, it has structural problems.** Groups are seeded greedily and
  depend on input order. Members can sit up to 2×tol apart (347 loci). Tolerant pairs get
  split (49 loci lost). Shifted members drag exact chains into the reference filter (45 lost).
  Chains shared exactly are labelled `tol10` (4,244 loci). 6,727 representatives have a splice
  site no other species supports.
- **gffcompare's own quirks.** With several inputs it silently drops same-species duplicate
  chains (`-D` is implied), and by default it merges exons separated by < 5 bp. Both need
  handling: restore duplicates by chain lookup, and pass `--no-merge`.

Decision: **exact SJs only**, with support defined by gffcompare. Exact matching removes every
tolerance problem by construction.

### Why "supported chains only" for grouping and linking (measured 2026-10-06)

Linking **all** models through any shared SJ creates giant components, even after the
reference filter: Earias 1,475 models / 230 species / 1,405 distinct chains / 287 kb in one
group; Lutraria 2,478 models / 5 species / 413 chains. Picking the representative by
"most SJs" over all models chooses a model with single-species junctions in 30–60 % of groups.
Restricting both grouping and representative choice to supported chains fixes both:

With the decided rules (exon-level same-strand reference filter per model, self excluded,
supported chains only):

| target | supported chains | genes | genes with 1 chain | max chains / gene | gene span p50 / max (bp) |
|---|---|---|---|---|---|
| Dermacentor andersoni | 116 | 107 | 101 | 4 | 3,378 / 151,026 |
| Drosophila melanogaster | 221 | 212 | 203 | 2 | 1,171 / 22,553 |
| Solanum lycopersicum | 170 | 160 | 152 | 3 | 1,397 / 85,031 |
| Earias clorana | 3,797 | 3,622 | 3,523 | 9 | 1,261 / 309,011 |
| Lutraria lutraria | 1,638 | 1,588 | 1,581 | 34 | 547 / 229,443 |

For comparison, `gffcompare_shared_chains.py` kept 115 chains on Dermacentor andersoni
(group-level filter). The 1 extra chain is rescued by the per-model filter: an overlapping
model is dropped, but the chain keeps ≥ 2 species.

### Why gffcompare does not do the grouping

gffcompare's tracking file groups **identical** intron chains only, and its XLOC loci group by
exon overlap. "Share ≥ 1 SJ" is neither, so the gene grouping is a union-find over an
SJ → chain index in Python. That is exact hashing with no thresholds. gffcompare remains the
authority on chain identity. Python re-derives the same partition and **exits 1 if they
disagree**, which guards against gffcompare version or option drift.

## New flow

```
PROJECTION.out.allmodels_raw  [meta(target_id, feature_type, decoy, gtype), allModels.<gtype>.raw.gff3]  × F·D·T
        ⋈ target_id
ch_targets                    [meta(id = target), fasta, gff3 | []]   (full samplesheet GFF3; [] for a pure target)
        │
        ▼
CURATE_MODELS (bin/curate_models.py)            one task per (target, feature_type, decoy)
  1. read models          species from `source=` (fallback: ID suffix after the last `|`); unique IDs asserted
  2. drop                 self (default) · single-exon · unstranded
  3. reference filter     drop models with an exon overlapping ANY ref exon entry, same strand; record the
                          overlapped ref types (if ref). Supported chains are counted before AND after this step
  4. exact chains         one GTF per species (integer IDs) → gffcompare --no-merge -i list → .tracking
                          cliques; verify vs Python chain keys; restore same-species duplicates
  5. supported chains     ≥ min_species distinct species
  6. genes                union-find over supported chains sharing ≥ 1 SJ
  7. representative       max chain species → max n_sj → max exonic length → min de → ID; ref_location (§2a)
  8. write                curated.gff3.gz · merged.gff3.gz (if ref) · report tables (observed / discarded by type / novel stats)
  9. validate             gffread on the curated GFF3 (exit ≠ 0 → task fails); merged: curated IDs must not
                          clash with reference IDs (exit 1), seqid blocks contiguous by construction
        │
        ▼
targets/<target>/curated/   <target>.curated.<gtype>.gff3.gz
                            <target>.curated.<gtype>.merged.gff3.gz     (annotated targets)
mqc_files → REPORTING (per-target report tables) and RUN_SUMMARY
```

## Design

### 1. `bin/curate_models.py` (new) — usable standalone and from Nextflow

```
curate_models.py --gff <fomo GFF3[.gz]> [--ref <reference GFF3/GTF[.gz]>] --target <species>
                 --prefix <out prefix> [--min-species 2] [--keep-self]
                 [--track lnc_RNA|mRNA|lnc_RNA.decoy|…] [--strand same|both]
                 [--gene-prefix <ID prefix>] [--gffcompare gffcompare] [--gffread gffread] [--skip-gffread]
                 [--keep-tmp] [--chains-tsv]

`--chains-tsv` writes an audit table, `<prefix>.chains.tsv.gz`: one row per supported chain
before the filter, with its species before and after, status (`kept` / `lost_ref`), overlapped
reference types, gene ID, representative, and every model. Off in the pipeline; used for the
phase-A cross-check.
```

- **Input.** Any fomo GFF3 with `gene` / `transcript` (or `lnc_RNA` / `mRNA`) / `exon` rows and
  `source=` on the transcript: `allModels.*.raw.gff3`, a `from_<source>` file, `top3`, or an old
  run whose IDs are not yet `|source`-qualified.
  - **Exit 1** on a duplicate transcript ID, an exon whose parent is missing, or exons on mixed
    seqids or strands within one transcript.
  - **An empty input is not an error.** It writes empty curated output, a merged file that equals
    the reference, and a funnel row of zeros (the "nothing projected" case, see CLAUDE.md
    merge/split contract).
- **Reference reader.** Every `exon` row of the reference, whatever its parent, from GFF3 or
  GTF, `.gz` OK.
  - Each exon carries a **reference type** for the report: the parent feature's type (`mRNA`,
    `lnc_RNA`, `antisense_RNA`, `transcript`, `tRNA`, `pseudogenic_transcript`, …). If that type
    is the generic `transcript` / `ncRNA` / `primary_transcript`, the type is refined by its
    biotype attribute (`transcript_biotype`, `biotype`, `ncrna_class`, `gbkey`), e.g.
    `transcript:misc_RNA`.
  - A transcript-like feature (it has a `Parent` gene or is itself the parent of CDS/UTR rows)
    with **no exon children** counts by its own span. CDS-only RefSeq models and a few
    `ncRNA_gene`-less records would otherwise be invisible.
  - Gene-level features (`gene`, `pseudogene`, `ncRNA_gene`, `region`) are not exons and are
    ignored. Their transcripts' exons carry them.
  - Overlap uses a sorted-starts + max-length bisect per `(chr, strand)`, exon vs exon, at least
    1 bp. With `--strand both`, the opposite strand is checked too.
- **gffcompare step.**
  - Write one GTF per species, with `transcript_id "t<index>"`. Never pass the real IDs: they
    contain `|`, which is the field separator inside `.tracking` `q` columns.
  - List the files in `inputs.txt` and run `gffcompare --no-merge -o cmp -i inputs.txt`.
  - Parse `cmp.tracking`: each line is one clique, at most one transcript per sample
    (`qN:g<i>|t<i>|…`, comma-separated if ever several).
  - **Checks (exit 1 on failure):**
    - each clique's members have one identical `(chr, strand, chain)`;
    - no chain appears in two cliques;
    - every multi-exon model's chain appears in some clique.
  - Then expand each clique to every model with that chain. This restores the duplicates `-D`
    dropped, so `supporting_models` is complete.
  - Validated on 82 targets: 0 check failures, 3,624 duplicates restored on one tick target
    alone; 524 inputs in < 1 min, ~1 GB RSS (Earias).
- **Genes.** Union-find keyed by SJ over supported chains only. Gene IDs are
  `<gene-prefix>_G<6 digits>`, numbered by `(seqid order of first appearance in input, start, strand)`
  so output is deterministic. The default `--gene-prefix` is `<prefix>`; Nextflow passes
  `FOMO_<target>_<gtype>` so IDs can never collide with reference IDs (asserted for the merged
  file).
- **Representative.** Key `(-n_sj, -n_species(chain, after the filter), -exonic_len, de, ID)`. A
  missing or NaN `de` sorts last. (Rev 3 = the phase-A key; rev 2 had species first.)

### 2. Output GFF3 (curated)

```
##gff-version 3
# curate_models: chains shared exactly by >=2 species (gffcompare --no-merge), linked by shared SJs; ref filter: same-strand overlap with any reference exon
<seq> fomo_curated gene        s e . ± . ID=FOMO_<t>_<gtype>_G000001;gene_biotype=fomo_lncRNA;n_chains=2;n_species_gene=3;ref_location=intergenic;curated_track=<gtype>
<seq> fomo_curated <tx type>   s e <score> ± . ID=<rep id>;Parent=FOMO_…_G000001;source=<sp>;…original attrs (de, AS, NM, …)…;
                                             transcript_biotype=fomo_lncRNA;n_supporting_species=3;supporting_species=A,B,C;supporting_models=id1|A,id2|B,id3|C;n_sj=4;
                                             ref_location=intergenic;ref_nearest=lnc_RNA;ref_distance=1010
<seq> fomo_curated exon        …           ID=<rep id>.exon1;Parent=<rep id>
```

### 2a. `ref_location` (rev 2)

Assigned to the representative. "Span" means a reference transcript's first-to-last exon (CDS/UTR
rows, or its own span, when it has no exons); it comes from the same spans index as the filter.
Those spans are built per (parent, seqid, strand), so reference rows with strand `.` / `?` count
on **both** strands. First match wins:

| class | rule |
|---|---|
| `antisense_exonic` | an exon overlaps a reference exon on the opposite strand |
| `intronic` | the whole span lies inside one same-strand reference transcript span (no exon overlap: guaranteed by the filter) |
| `sense_span_overlap` | overlaps a same-strand reference span without being contained (straddles a span end, or contains a reference transcript inside its intron) |
| `antisense_intronic` | overlaps an opposite-strand reference span, no exon overlap |
| `divergent` | no overlap; head-to-head with an opposite-strand reference transcript whose 5' end is ≤ `--divergent-dist` (default 1000 bp) upstream of ours |
| `intergenic` | none of the above (tail-to-tail neighbours are intergenic) |
| `no_reference` | run without `--ref` |

- `ref_nearest`: the type label of the nearest reference transcript span, either strand.
- `ref_distance`: in bp, 0 when overlapping. Both are `.` when the seqid has no reference
  transcript, and absent without `--ref`.
- Cost: only representatives are queried, and the span index stores type labels only (no
  transcript IDs). Mercenaria still peaks at 772 MB.
- With `--strand both` the filter removes opposite-strand exon overlaps too, so
  `antisense_exonic` cannot occur.

- **Escaping.** Attribute values are GFF3-escaped (`; = & ,` → `%XX`) **inside** each list
  element; the commas between elements stay literal.
- **Sorting and layout.** Records are grouped by seqid, then sorted by start. Gene → transcript →
  exons are contiguous.
- **Transcript type.** Kept from the input row (`transcript` for current fomo output).
  `gene_biotype` and `transcript_biotype` come from the track: `lnc_RNA` → `fomo_lncRNA`,
  `mRNA` → `fomo_mRNA`, decoy → `fomo_decoy_lncRNA` / `fomo_decoy_mRNA`. The script takes the
  track from `--track` (Nextflow passes `meta.gtype`); standalone it guesses from the file name
  (`.lnc_RNA.` / `.mRNA.` / `.decoy.`), falling back to `fomo_lncRNA`.
- **Compression.** Written with Python `gzip` (plain gzip, `.gff3.gz`). Every GFF3 consumer in
  the pipeline reads `.gz` (CLAUDE.md, plans/17).

### 3. Merged GFF3 (curated + reference)

- **Content.** The reference is copied **verbatim**: every directive, comment and feature line,
  all attributes kept, decompressed. Each curated gene block is inserted **after the last
  reference line of its seqid**; seqids absent from the reference are appended at the end. This
  keeps every seqid contiguous, which `gff-feature-stats` requires above its 200 k-line batch
  threshold (CLAUDE.md). It avoids re-sorting reference hierarchies (RefSeq children are not
  always adjacent to parents) and avoids gffread rewriting the reference: gffread drops
  `region` rows and renames attributes unless told otherwise.
- **Not re-validated with gffread.** The reference part is verbatim, and gffread hard-errors on
  some real references: RefSeq plant mitochondrial trans-spliced models with strand `?`, seen on
  10 solanaceae/apiaceae targets in phase A. The curated part is validated with gffread
  on its own. Curated IDs are checked against every reference `ID=` in the same streaming pass
  that finds each seqid's last line.
- **No same-strand overlap with any reference exon, by construction.** Curated models may
  still sit inside reference introns or be antisense to a reference exon; that is the chosen
  rule.
- **A pure target (no reference)** gets no merged file. `curated.gff3.gz` is its deliverable.
- **A GTF reference** (standalone use only; the samplesheet requires GFF3) is converted with
  `gffread --keep-genes -F` before merging, and the script logs that it did so.

### 4. Report tables (one set per track)

The script writes `<prefix>.curation.tsv`: long format, `table	metric	value`, one row per
number, so MultiQC adapters can reshape it without re-parsing GFF3. It also writes
`<prefix>.curation.json` with the same numbers. Three tables:

**(1) Observed** — multi-species models before the reference filter:
- input: `n_models`, `self_excluded`, `single_exon`, `multi_exon_models`;
- models in supported chains: `models_in_supported_chains`;
- `chains` (distinct), `supported_chains` (≥ min species);
- `genes` (supported chains linked by SJs);
- `supported_chains_by_n_species` (2 / 3 / 4 / ≥ 5);
- `dup_restored` (gffcompare `-D` duplicates restored).

**(2) Discarded by reference overlap** (only with a reference):
- `models_removed`;
- `supported_chains_lost`: supported before the filter, < min species after;
- `supported_chains_kept_after_removal`: a model was removed but support stays ≥ min;
- `genes_lost`.

Then the same counts **per overlapped reference type** (`mRNA`, `lnc_RNA`, `antisense_RNA`,
`transcript:misc_RNA`, `tRNA`, …):
- a model or chain overlapping several types counts once under each (`by_type`);
- it also counts once under a single primary type (`by_primary_type`), chosen by priority
  `mRNA > lnc_RNA > other lnc-like (antisense_RNA, lncRNA biotypes) > transcript:* > small RNAs >
  pseudogenic > other`;
- the primary counts add up to the total.

**(3) Novel models** — the kept genes / representatives:
- `genes`, `genes_by_n_species` (2 / 3 / 4 / ≥ 5), `genes_by_n_chains` (1 / 2 / ≥ 3);
- representative `n_exons` (min / median / mean / max);
- `transcript_length` (exonic bp; min / p25 / median / mean / p75 / max);
- `gene_span` (genomic bp, same stats);
- `intron_length` (same stats over all representative introns);
- `de` (median / mean);
- `genes_location:<class>` for all 7 classes of §2a, every class always written (0 when empty).
  They sum to `genes`. Rev 2 replaced the phase-A `genes_antisense_to_ref_exon` /
  `genes_within_ref_transcript_span_same_strand`.

The pipeline (phase B) turns (1)–(3) into three MultiQC custom-content tables per target. The
run summary (phase C) takes per-target totals, plus the decoy-normalised FDR.

### 5. Nextflow wiring

- **`modules/local/curate_models.nf`** (new): `CURATE_MODELS`, label `process_low` (sized
  after the real-data dry run below).
  - Input `tuple val(meta), path(gff3), path(ref)`, where `ref` may be `[]`.
  - Outputs:
    - `curated` (`*.curated.*.gff3.gz`, always);
    - `merged` (`*.merged.gff3.gz`, `optional: true`);
    - `report` (`*.curation.{tsv,json}`) → `*_TO_MQC` adapter (phase B);
    - versions topic: python, gffcompare, gffread.
  - Container: one **Wave** image built from conda `python=3.11`, `bioconda::gffcompare=0.12.6`,
    `bioconda::gffread=0.12.7`. Those are the versions already pinned in
    `modules/nf-core/{gffcompare,gffread}`, so curation and benchmarking agree. The 82-target
    validation used 0.12.10, so re-run the CLI cross-check under 0.12.6. Add
    `environment.yml`. Use the Seqera MCP tools per CLAUDE.md.
- **`subworkflows/local/curation.nf`** (new): a **top-level** `CURATION` called from
  `workflows/fomo.nf`, not nested in PROJECTION (CLAUDE.md multi-target rule 1).
  - `take: allmodels_raw, ch_targets`.
  - Join: `allmodels_raw.map { m, g -> [m.target_id, m, g] }.combine(ch_targets.map { m, _f, g -> [m.id, g ?: []] }, by: 0)`.
    Every target is exactly one `ch_targets` row, so the combine is a keyed 1:1. It covers both
    pure targets (`[]`) and `both` targets (full samplesheet GFF3, **not** the filtered ★
    annotation, which is lncRNA-only and longest-isoform).
  - Emit `curated`, `merged`, `mqc_files` (meta carries `target_id` → per-target report).
- **`workflows/fomo.nf`**: call `CURATION` when `params.curate`, and mix its `mqc_files` into
  `ch_all_mqc`.
- **`conf/modules.config`**:
  - `ext.prefix = { "${meta.target_id}.curated.${meta.gtype}" }`;
  - `ext.args = { "--target ${meta.target_id} --min-species ${params.curate_min_species} --gene-prefix FOMO_${meta.target_id}_${meta.gtype}" }`;
  - publishDir `targets/<target>/curated/`, anchored `'.*:CURATE_MODELS$'`, not publishing
    `_mqc.tsv` or `versions.yml`.
- **Params**: `curate = true`, `curate_min_species = 2`. Add both to **`nextflow_schema.json`**
  (`failUnrecognisedParams`). `-profile test` keeps `curate = true`.
- **MultiQC**:
  - `assets/multiqc/sections.yml`: a `curation` custom_data table, `sp: fn: '*_curation_mqc.tsv'`;
  - `main.yml`: a `report_section_order` slot after the projection sections, and
    `_curation_mqc` added to `extra_fn_clean_exts`;
  - no `# id:` header in the TSV (rule 1 of "Adding a new custom-content section").
- **Run summary** (phase C): add `curated_genes` (lncRNA) and, with decoys,
  `curated_decoy_genes` / `curated_fdr = decoy / real` per target to `run_summary_tables.py`.
- **Docs**: CLAUDE.md gets the `targets/<target>/curated/` tree, a `CURATION` row in the
  subworkflow map, and `F·D·T` `CURATE_MODELS` in the task-count paragraph. Add README outputs.

### 6. Consequences to know before running

- **mRNA track.** On an annotated target, nearly every projected mRNA lands on an orthologous
  reference mRNA, so the any-exon filter removes it. The curated mRNA set is then *novel* coding
  loci only. On a pure target it is every supported coding chain.
- **Decoy track.** Each source's decoys are relocated independently into its *own* intergenic
  space. Two sources rarely project a decoy onto the same exact target chain, so curated decoys
  should be near 0. That is the point: it bounds the false-positive rate of the curated set.
  With `decoy_cap = 1000` per source, decoy and real counts are on different scales, so the FDR
  in the run summary is normalised per input model.
- **Disk.** Each merged file contains the full reference, once per track (up to 4 per target
  with both flags), hence gzip.
- **Counts differ from `curate_shared_chains.py`.** That script used tol 10, a group-level
  reference filter, one model per tolerance group and the exon filter on all members. The
  cross-check below must explain every difference.

## Implementation phases

- **A — `bin/curate_models.py` + fixtures** (no Nextflow).
  - Fixtures in `tests/curate_models/`: small GFF3s plus expected outputs, run by `run.sh`.
  - Cases:
    1. chain shared by 2 species → kept;
    2. same chain twice in one species → not supported;
    3. two supported chains sharing 1 SJ → one gene, rep = more SJs;
    4. two supported chains linked only through an unsupported model → two genes;
    5. same-strand exon overlap removes one model and drops support below 2 → chain gone;
    6. opposite-strand overlap → kept;
    7. intronic (no exon overlap) → kept;
    8. representative ties at each key (n_sj, species, length, de, ID);
    9. IDs containing `|`, `,`, `;`;
    10. single-exon only;
    11. empty input;
    12. no reference;
    13. GTF reference;
    14. self models excluded;
    15. merged: unique IDs, seqid contiguity, gffread passes;
    16. reference types: exon of tRNA / antisense_RNA / misc_RNA `transcript` / exon-less
        CDS-only model each remove a model and are reported under their type; a model
        overlapping two types counts once per type and once under its primary type.
- **B — Nextflow module + subworkflow + params + publish + MultiQC funnel.**
- **C — Run summary columns (curated genes, decoy FDR) + docs.**

## Files to touch

| File | Change |
|---|---|
| `bin/curate_models.py` | new |
| `tests/curate_models/` | new fixtures + `run.sh` |
| `modules/local/curate_models.nf`, `modules/local/curate_models/environment.yml` (or alongside) | new |
| `subworkflows/local/curation.nf` | new |
| `workflows/fomo.nf` | call `CURATION`, mix `mqc_files` |
| `conf/modules.config` | ext.prefix / ext.args / publishDir for `CURATE_MODELS` |
| `nextflow.config`, `nextflow_schema.json` | `curate`, `curate_min_species` |
| `assets/multiqc/sections.yml`, `assets/multiqc/main.yml` | curation funnel section |
| `bin/run_summary_tables.py`, `assets/multiqc/run_summary_sections.yml` | phase C |
| `CLAUDE.md`, `README.md` | layout, subworkflow map, task counts |

## Verification

1. **Fixtures**: `tests/curate_models/run.sh` passes. Reruns are byte-identical.
2. **Cross-check vs the 2026-10-05 validation.**
   - Run the CLI on the same 82 targets (old-format IDs; `--ref` = the reference from each
     `curated_set/curate_commands.txt`).
   - (a) The supported-chain set must equal `gffcompare_shared_chains.py`'s kept set, plus
     chains rescued by the model-level filter. List every extra chain with the model removed.
   - (b) Every exact-labelled locus of the production `shared_chain.gff3` must be inside a
     curated gene, unless the new filter explains its absence.
   - (c) Report gene counts, the chains-per-gene distribution and the max span next to the
     table above.
3. **gffcompare version**: repeat 2(a) on 5 targets inside the Wave container (0.12.6).
   Results must be identical.
4. **`-profile test,singularity`**, default and with `--include_mrna --include_decoy`:
   - one curated GFF3 per (target, gtype), and a merged GFF3 for each `both` target;
   - `gffread` and `gff-feature-stats` accept every output;
   - the funnel section renders once per report;
   - no `targets/null/`.
5. **Stub run** and `nextflow lint .`: no new errors beyond the 7 known on `main`.
6. **Real data**: one acari and one mollusca target through the module (`-resume` on an
   existing work dir, or the CLI in the container) for runtime/RSS, used to size the label.

## Resolved questions (2026-10-06, second round)

1. **Reference types for the filter**: **any exon entry**, every biotype, antisense_RNA
   included. Same strand, as decided in the first round. `--strand both` is available if this
   should become "any strand" (to confirm).
2. **Compression**: outputs are **gzipped** (`.gff3.gz`).
3. **Biotype labels**: `fomo_`-prefixed: `fomo_lncRNA`, `fomo_mRNA`, `fomo_decoy_lncRNA`,
   `fomo_decoy_mRNA`.
4. **Self support**: never counted. Expected to be rare after the filter.
5. **Report**: three tables per track. (1) Multi-species models observed, (2) discarded by
   reference overlap with the overlapped types, (3) summary stats of the novel models. See §4.

## Out of scope / follow-ups

- Curating `top3`.
- Tolerant SJ matching. Rejected in the 2026-10-05 review: greedy seeding, mislabelling,
  2×tol drift.
- Benchmarking the curated set against the reference (gffcompare Sn/Pr of supported chains
  *before* the reference filter would measure how often cross-species support recovers known
  lncRNAs). Natural next step.
- Retiring `tools/fomo-curate/curate_shared_chains.py`. Keep it until phase C is verified;
  `gffcompare_shared_chains.py` stays as the reference for the cross-check.

## Phase A results (2026-10-06)

The CLI was run on the 82 targets of the 2026-10-05 review: all of acari (14), apiaceae (5),
solanaceae (13) and drosophilidae (44), plus 3 lepidoptera and 3 mollusca. Input was the existing
`allModels.lnc_RNA.raw.gff3` (old-format IDs, without `|source`; the new format is covered by
fixtures) and the reference from each `curated_set/curate_commands.txt`. **0 errors**.

**Correctness checks, all 82 targets:**
- Supported chains before the filter equal an independent Python hashing of the raw models,
  and each chain lists the same models.
- Every gene's representative is consistent with the rule: max SJs, then species; its
  `supporting_models` are exactly the surviving models of its chain.
- Every production `match=exact` locus is a supported chain (0 missing).
- Old-rule vs new-rule differences are all explained (below). 0 unexplained.
- Byte-identical results before and after the memory rewrite.

| | acari | apiaceae | solanaceae | drosophilidae | lepidoptera (3) | mollusca (3) | all |
|---|---|---|---|---|---|---|---|
| multi-exon non-self models | 141,754 | 19,130 | 250,761 | 322,907 | 243,952 | 166,116 | 1,144,620 |
| (1) supported chains observed | 2,190 | 431 | 10,155 | 16,242 | 16,853 | 5,056 | 50,927 |
| (2) supported chains lost to the reference | 833 | 347 | 6,198 | 5,066 | 5,499 | 1,803 | 19,746 |
| (3) novel genes | 1,307 | 84 | 3,623 | 10,658 | 10,948 | 3,178 | 29,798 |
| — 2 / 3 / 4 / ≥ 5 species | 1239/60/6/2 | 67/14/3/0 | 3231/218/100/74 | 8234/2004/299/121 | 7191/1900/753/1104 | 2281/520/206/171 | 22243/4716/1367/1472 |
| — 1 / 2 / ≥ 3 chains | 1275/20/12 | 84/0/0 | 3395/170/58 | 10216/393/49 | 10706/171/71 | 3160/9/9 | 28836/763/199 |
| — antisense to a ref exon | 645 | 26 | 678 | 3,005 | 5,064 | 320 | 9,738 |
| — within a same-strand ref transcript span | 170 | 1 | 378 | 669 | 2,191 | 946 | 4,355 |

Chains lost to the reference, by primary type:

| primary type | chains lost |
|---|---|
| mRNA | 10,447 |
| lnc_RNA | 8,706 |
| pseudogene | 231 |
| antisense_RNA | 106 |
| transcript:misc_RNA | 90 |
| gene:lncRNA (exons parented by the gene directly) | 58 |
| gene:protein_coding | 32 |
| rRNA, snRNA, snoRNA, tRNA, miRNA | 26, 13, 10, 8, 8 |

Median novel representative (median of the per-target medians): 2 exons, 547 nt exonic,
1.15 kb span, 247-nt introns, de 0.138.

**Against the old exact rule** (`gffcompare_shared_chains.py`: group-level, mRNA/lncRNA exons
only):
- 31,566 → 31,181 kept chains.
- **503 old-only** chains: all hit a reference type the old reader ignored (pseudogene 231,
  antisense_RNA 106, misc_RNA `transcript` 89, rRNA 26, …).
- **118 new-only** chains: all rescued by the per-model filter (an overlapping model dropped,
  ≥ 2 species left).

**Against the production tol10 output:**
- Of 25,679 `exact` loci, 25,221 are kept and 458 are lost to the new reference types.
- Of 11,265 `tol10` loci, 4,208 have a representative chain that is exactly supported and kept.
  The rest have no exact cross-species chain and are not reportable under this rule, as decided.
- Representatives can differ from production: the order is now SJs → species → **length** → de.
  Production picked by splice-site votes, then de.

**Resources:** median 4 s and max 28 s per target. Peak RSS max 772 MB (Mercenaria: 191 k
models, 176 species, 1.14 M reference exons). `process_low` (4 GB) is enough.

Making that fit took three changes:
- Models are stored compactly (slots, interned strings, raw line parsed only for
  representatives): 1.08 GB → 0.31 GB.
- The reference index uses typed arrays with labels resolved while reading: 1.2 GB → 0.5 GB.
- The merged file is no longer re-parsed by gffread.

**Not yet exercised on real data:**
- mRNA and decoy tracks (labels covered by fixtures);
- (new-format `|source` IDs turned out to be present already in the drosophilidae and mollusca
  runs, e.g. `lnc_RNA1971|Drosophila_carrolli_2789720`, so they were exercised);
- gffcompare 0.12.6 (the pipeline pin; the runs used 0.12.10). Phase B verification step 3.

### Revision 2 (rep by species, ref_location) — 2026-10-06

**Changes:**
- Representative = the chain shared by the most species first, then SJs, length, de, ID.
- `ref_location` / `ref_nearest` / `ref_distance` (§2a).
- `genes_location:*` in the `novel` table.
- New CLI option `--divergent-dist`.

**Fixtures: 61 checks pass** (43 before).
- New case 17: a 1-SJ chain in 3 species beats a 3-SJ chain in 2 species that shares an SJ with
  it.
- New case 18: each location class on a new reference `chr5`. It includes divergent on both
  strands, `+` 200 bp → divergent vs 1,010 bp → intergenic, and tail-to-tail → intergenic. It
  checks `ref_distance` / `ref_nearest`, that gene and transcript carry the same class, the GTF
  reference, and `no_reference`.
- No earlier expectation changed. In case 3 both chains have 2 species, so SJs still decide; in
  case 8 R2 already won on species.

**10 random targets** (`random.seed(21)` over the 82), plus Mercenaria for the RSS re-check.
Checks a/b/c/d (d now species-first) all pass, the location counts sum to `genes`, and the gene
sets (by chain composition) are identical to phase A on every target:

| target | genes | reps changed | antisense_exonic | intronic | sense_span_overlap | antisense_intronic | divergent | intergenic | wall s | RSS MB |
|---|---|---|---|---|---|---|---|---|---|---|
| Lycium ferocissimum | 216 | 3 | 45 | 9 | 9 | 8 | 6 | 139 | 8.8 | 179 |
| Drosophila mojavensis | 101 | 0 | 35 | 3 | 3 | 14 | 5 | 41 | 3.6 | 86 |
| Lutraria lutraria | 1,582 | 6 | 189 | 116 | 71 | 103 | 1 | 1,102 | 8.1 | 196 |
| Drosophila bipectinata | 250 | 1 | 69 | 7 | 4 | 43 | 14 | 113 | 3.9 | 92 |
| Drosophila rhopaloa | 495 | 5 | 137 | 15 | 5 | 89 | 18 | 231 | 3.9 | 91 |
| Solanum lycopersicum | 147 | 0 | 48 | 4 | 10 | 9 | 2 | 74 | 8.8 | 177 |
| Drosophila pseudoobscura | 200 | 2 | 50 | 11 | 3 | 38 | 2 | 96 | 4.6 | 104 |
| Drosophila subobscura | 226 | 2 | 74 | 4 | 3 | 37 | 9 | 99 | 3.7 | 84 |
| Nicotiana sylvestris | 425 | 0 | 41 | 7 | 29 | 33 | 8 | 307 | 10.4 | 203 |
| Drosophila simulans | 354 | 1 | 110 | 16 | 2 | 55 | 15 | 156 | 4.5 | 100 |
| **10 targets** | **3,996** | **20 (0.5 %)** | 798 (20 %) | 192 (5 %) | 139 (3 %) | 429 (11 %) | 80 (2 %) | 2,358 (59 %) | | |
| Mercenaria mercenaria (RSS check) | 1,474 | 4 | 103 | 549 | 108 | 370 | 15 | 329 | 27.0 | 772 |

Every changed representative is in a multi-chain gene: a 2-SJ chain in 2 species is replaced by a
1-SJ chain in 3–7 species. Examples (old → new rep, species / SJs):
- Lutraria: `ENSQQCT00005065900|Batillaria_attramentaria_370345` (2 sp / 2 SJ) →
  `ENSTEST00000077716|Ostrea_denselamellosa_74434` (7 sp / 1 SJ).
- D. pseudoobscura: `rna-XR_001768618.3|Drosophila_biarmipes_125945` (2 / 2) →
  `rna-XR_001772318.2|Drosophila_eugracilis_29029` (7 / 1).
- D. bipectinata: `lnc_RNA1369|Drosophila_carrolli_2789720` (2 / 2) →
  `rna-NR_124784.1|Drosophila_melanogaster_7227` (3 / 1).

Phase-A `genes_antisense_to_ref_exon` (33 % over 82 targets) corresponds to `antisense_exonic`
here. The phase-A "within a same-strand span" count is now split into `intronic` (contained)
and `sense_span_overlap` (partial).

### Revision 3 (2026-10-06): representative back to SJs first

At the user's request the representative key is again `(-n_sj, -n_species, -exonic_len, de, ID)`,
so a gene reports its most complete supported chain. Rev 2 had species first: on 10 targets it
changed 20 / 3,996 representatives (all multi-chain genes; a 2-SJ chain in 2 species replaced by a
1-SJ chain in 3–7 species). Those revert to the phase-A choice. Gene grouping and `ref_location`
(rev 2) are unchanged.

Fixture 17 is inverted: the 3-SJ chain in 2 species beats the 1-SJ chain in 3 species (gene
`n_chains=2`, `n_species_gene=5`).
