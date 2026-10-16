# 25 — Curated GFF3 format: content-hashed IDs, one transcript per supported chain

> **Status:** Implemented · **Date:** 2026-10-08
> **Decision:** Curated IDs are content hashes of the gene's intron chains; a gene has one transcript per supported chain, ranked by support.

## Status: implemented (2026-10-08), not yet committed — see **Results**

## Decisions (2026-10-08)

| Question | Decision |
|---|---|
| Transcripts per gene | **One per supported intron chain** (was: one representative per gene). A gene with `n_chains=2` gets `t1` and `t2`. |
| Representative of a chain | The model with the **longest spliced (exonic) length** among the chain's supporting models (after the reference filter) → lowest `de` → model ID. Unchanged in effect: within one chain the current key already reduces to this. |
| Transcript order (`t1`, `t2`, …) | **By support**: most supporting species → most SJs → longest spliced length → leftmost start → rightmost end → representative model ID. `t1` is the best-supported isoform, not necessarily the longest. |
| ID scheme | Content hash of the **gene's intron chains on the target** (not of source transcript IDs). One hash `H` per gene, shared by its gene, transcripts and exons. |
| Gene ID | `fomo_<tag>_g_<H>` |
| Transcript ID | `fomo_<tag>_t<n>_<H>` |
| Exon ID | `fomo_<tag>_t<n>e<k>_<H>`, `k` numbered **5′→3′** (on `-` the rightmost exon is `e1`) |
| Column 3 (transcript) | The pipeline feature type verbatim: `lnc_RNA` or `mRNA` (SO terms). **Decoys use the same type** as their real track. |
| Column 6 | `.` everywhere (was the minimap2 MAPQ on the transcript) |
| Alignment metrics | **Dropped** from the curated GFF3 (`query_length`, `aligned_bases`, `pct_*`, `projected_length`, `NM`, `AS`, `de`). They stay in `annotation/<target>.from_<source>.<gtype>.projected.gff3`, keyed by the `tid\|species` IDs listed in `supporting_transcripts`. |
| `ref_*` attributes | **Gene only.** Missing values are **omitted**, never written as `.`. |
| Biotype | `biotype=` on gene **and** transcript: `lncRNA`, `mRNA`, `decoy_lncRNA`, `decoy_mRNA` (was `fomo_lncRNA` etc.; provenance is column 2, `fomo`). |

## Why

The current curated GFF3 (example: `results_OctVul`, `G000007`) has five problems:

- **Transcript IDs are another species' ID.** Example: `ID=ENSMQKT00000025394|Eumandya_parva_2302672`.
- **Gene IDs are ordinals.** `FOMO_<target>_lnc_RNA_G000007` shifts whenever anything upstream
  changes, so IDs cannot be compared across runs.
- **The transcript line is overloaded.** It has about 30 attributes: alignment metrics that
  describe only one of the supporting alignments, three redundant biotype keys, the `ref_*` tags
  duplicated from the gene, and `source=`.
- **Supported isoforms are dropped.** A gene whose chains are all supported by ≥ 2 species
  reports only one of them, so the other supported isoforms are discarded. This happens for real:
  the apiaceae run (`results_new`) has up to 3 two-chain genes per target.
- **Column 6 holds the minimap2 MAPQ (60)**, which reads as a model score.

**Why hash the chains and not the source transcript IDs.** A curated model *is* its exact intron
chain, so the chain is what the ID should identify. The supporting-ID set changes under three
conditions that do not change the model:

- a source ID's version bumps (NCBI `XR_….1` → `.2`);
- a source species is added to or removed from the run;
- the reference annotation changes, which changes which models survive the overlap filter.

A chain hash is stable as long as the target assembly is. So the same model gets the same ID in
every run, whatever the source set, and runs can be compared by ID.

**Why one hash for the whole gene, with ordinals.** It groups a gene's rows visually and under
grep. If a new chain joins a gene, `t1`/`t2` can reorder, but `H` changes too, because the chain
set changed, so the ordinals cost no stability that the gene hash had not already lost.

## Design

All changes are in `bin/curate_models.py`, except where noted under **Files**.

### 1. Hash

```
track_key = the --track string (lnc_RNA | mRNA | lnc_RNA.decoy | mRNA.decoy)
chain_str = ",".join(f"{a}-{b}" for a, b in introns)          # genomic order, 1-based intron coords
gene_key  = f"{track_key}|{seqid}|{strand}|" + ";".join(sorted(chain_str for each chain of the gene))
H         = sha1(gene_key).hexdigest()[:12]
```

- **Track in the key.** A chain curated in both the lncRNA and the mRNA track must not collide if
  the two files are ever concatenated. The `<tag>` prefix separates them too, but keeping the
  track in the key as well costs nothing.
- **Inputs.** The target id is **not** in the key: curated files of different targets live on
  different genomes and are never merged.
- **Collisions.** 12 hex characters is 48 bits, so a collision among 10⁵ genes has probability
  ~2·10⁻⁵. The script **exits 1** on a duplicate `H` within a file (naming both genes); it never
  silently suffixes.

### 2. `<tag>` and the prefix option

`<tag>` is derived from the track: `lncRNA`, `mRNA`, `decoy_lncRNA` or `decoy_mRNA`. It is the
same token as `biotype=`. `track_info()` therefore returns it without the `fomo_` prefix, and the
ID prefix is `fomo_<tag>`.

`--gene-prefix` is **replaced** by `--id-prefix` (default `fomo_<tag>`). The option is kept only
as the escape hatch for the merged-file ID clash check (`write_merged` exits 1 if a curated ID
already exists in the reference). In `conf/modules.config`, `CURATE_MODELS`' `ext.args` drops
`--gene-prefix FOMO_${target}_…`. The default is correct, and the target id no longer belongs in
the ID.

### 3. Genes and transcripts (`main`, after the reference filter)

For each gene (connected component of supported chains, `link_genes` unchanged):

```
for each chain k in the gene:
    cand_k = models with chain k that survived the reference filter
    rep_k  = min(cand_k, key=(-exonic_len, de, id))
    sup_k  = sorted(cand_k, key=(species, id))           # supporting_transcripts
    sp_k   = sp_post[k]                                  # supporting_species
transcripts = sorted(chains, key=(-len(sp_k), -len(k), -exonic_len(rep_k),
                                  rep_k.start, rep_k.end, rep_k.id))      → t1, t2, …
```

The gene **span** is the union of its transcripts' spans. Genes are sorted by
`(seqid order, gene start, strand, gene end, H)`.

### 4. Gene location against the reference

`ref_location()` runs **once per gene**, on a pseudo-model: the gene's seqid and strand plus the
**merged union of all its transcripts' exons** (overlapping or adjacent exons merged into a
sorted, disjoint list, so `ex[0][0]` / `ex[-1][1]` are the gene span).

- It is the gene's footprint, so the class answers "where is this gene", which is what a gene
  attribute should mean.
- For a single-transcript gene the footprint is exactly today's representative, so **every
  `n_chains=1` gene keeps its current `ref_location` and `ref_*` tags**. This makes it a
  regression check (see Verification).
- `check_location()` receives the gene ID for its error message.

### 5. Output rows

```
gene      fomo  gene        S  E  .  ±  .  ID=fomo_<tag>_g_H;biotype=<tag>;n_chains=N;n_species=M;ref_location=…[;ref_gene_id=…][;ref_gene_name=…][;ref_gene_biotype=…][;ref_gene_orientation=…][;ref_distance=…]
per t<n>  fomo  <lnc_RNA|mRNA> s  e  .  ±  .  ID=fomo_<tag>_t<n>_H;Parent=fomo_<tag>_g_H;biotype=<tag>;representative_transcript=<tid|species>;supporting_transcripts=<tid|species>,…;supporting_species=<sp>,…;n_sj=J
per exon  fomo  exon        a  b  .  ±  .  ID=fomo_<tag>_t<n>e<k>_H;Parent=fomo_<tag>_t<n>_H
```

- **Row order.** Gene, then `t1` and its exons, then `t2` and its exons, and so on.
- **Exon rows.** Written in genomic order; only the exon *number* `k` follows 5′→3′. A
  `-`-strand transcript therefore lists `t1e3, t1e2, t1e1`, as Ensembl GFF3 does.
- **One exon row per transcript**, even when `t1` and `t2` share an exon. No multi-`Parent`
  exons.
- **`n_species` (gene)** is the union of species over **all** the gene's chains.
  **`supporting_species` (transcript)** is that transcript's chain only. They are equal when
  `n_chains=1`.
- **`ref_*` keys.**
  - `no_reference`: only `ref_location=no_reference`.
  - Seqid with no reference transcript (`ref is None`): `ref_location=intergenic` and **no**
    `ref_*` keys.
  - `ref_gene_name` is omitted when the reference gene has no name (today `.`).
- **Escaping.** Values go through `esc()` as today, so commas inside source IDs stay `%2C` in
  `supporting_transcripts`.
- **Header comment.** The `# curate_models:` line changes "one representative per gene" to "one
  transcript per supported chain (ordered by support)", and reports genes **and** transcripts.

### 6. Report (`curation.{tsv,json}`) and MultiQC

- **`novel`.**
  - Add `transcripts` and `transcripts_per_gene_mean`.
  - `genes_<n>sp` counts genes by **gene-level** `n_species` (union). Today it counts the
    representative chain's species; the two differ only for multi-chain genes.
  - `n_exons`, `transcript_length`, `intron_length`, `supporting_models` and `de` medians are over
    **all transcripts**. `gene_span` is over genes (union span).
- **`bin/curation_to_mqc.py`.** Add `transcripts` to `NOVEL`.
- **`assets/multiqc/sections.yml`.** Rewrite the `curation_novel` description: "represented by
  its chain with the most splice junctions" becomes "one transcript per supported chain".
- **`bin/run_summary_tables.py`.** No change: it reads `genes`, whose meaning is unchanged.

### 7. `--chains-tsv`

A `transcript_id` column (`fomo_<tag>_t<n>_H` for every kept chain; `.` for `lost_ref`) is added
before `representative`, which is kept: it names the model whose exons the transcript carries.
`gene_id` takes the new gene ID.

### 8. Merged file

`new_ids` covers every gene, transcript and exon ID. Nothing else changes in `write_merged`.

## Example

`results_OctVul` `G000007` (one chain), before and after:

```
OX597815.1  fomo  gene     51846883  51848249  .  +  .  ID=fomo_lncRNA_g_<H>;biotype=lncRNA;n_chains=1;n_species=2;ref_location=antisense_exonic;ref_gene_id=gene:ENSVUCG00010032847;ref_gene_name=MYEF2;ref_gene_biotype=protein_coding;ref_gene_orientation=antisense;ref_distance=0
OX597815.1  fomo  lnc_RNA  51846883  51848249  .  +  .  ID=fomo_lncRNA_t1_<H>;Parent=fomo_lncRNA_g_<H>;biotype=lncRNA;representative_transcript=ENSMQKT00000025394|Eumandya_parva_2302672;supporting_transcripts=ENSMQKT00000025394|Eumandya_parva_2302672,ENSDVDT00000015758|Sepia_lycidas_296134;supporting_species=Eumandya_parva_2302672,Sepia_lycidas_296134;n_sj=1
OX597815.1  fomo  exon     51846883  51847268  .  +  .  ID=fomo_lncRNA_t1e1_<H>;Parent=fomo_lncRNA_t1_<H>
OX597815.1  fomo  exon     51848149  51848249  .  +  .  ID=fomo_lncRNA_t1e2_<H>;Parent=fomo_lncRNA_t1_<H>
```

## Files

| File | Change |
|---|---|
| `bin/curate_models.py` | §1–8; module docstring (Method 5–7, Output) |
| `bin/curation_to_mqc.py` | `transcripts` in `NOVEL` |
| `assets/multiqc/sections.yml` | `curation_novel` description |
| `conf/modules.config` | `CURATE_MODELS` `ext.args`: drop `--gene-prefix …` |
| `tests/curate_models/test_curate_models.py` | update and extend (below) |
| `CLAUDE.md`, `README.md` | the CURATION row / curation step: "one representative per gene" → "one transcript per supported chain"; the ID scheme |

## Tests (`tests/curate_models/test_curate_models.py`)

**Update existing cases.**

- **Case 3** (gene with 2 chains). Now **two** transcripts, ordered by support. The gene has
  `n_chains=2;n_species=4`, and each transcript lists its own chain's
  `supporting_transcripts` / `supporting_species`.
- **Case 17.** Same, under the new attribute names.
- **Case 9.** The escaped comma survives into `supporting_transcripts`.
- **Biotype checks.** `biotype=lncRNA` / `decoy_lncRNA` / `mRNA`, and column 3 is
  `lnc_RNA` / `lnc_RNA` / `mRNA`.
- **Clash test.** The reference carries `ID=fomo_lncRNA_g_<H>` of a real curated gene: exit 1.
  With `--id-prefix X`: passes.
- **Tag checks** (`TAG`, 19): read `ref_*` from the **gene** row.

**New cases.**

- **Stability.** The same models with the input file shuffled → byte-identical curated GFF3.
  Adding an unrelated source's models elsewhere → the IDs of untouched genes are unchanged.
- **Ordering.**
  - Two chains, where the one with more species has fewer SJs: it is `t1`.
  - Equal species: more SJs first.
  - Then longer first.
- **Representative.** Within a chain, the longest spliced model is chosen even when another
  model has a lower `de`.
- **Exon IDs.** Unique across the file. On a `-` strand transcript, `e1` is the rightmost exon,
  and rows are in genomic order.
- **Every ID unique** across the file (gene ∪ transcript ∪ exon).
- **Missing reference keys.**
  - Seqid without reference transcripts: `ref_location=intergenic`, no `ref_*` keys.
  - Unnamed reference gene: no `ref_gene_name` key.
  - `no_reference`: only `ref_location`.
- **Multi-chain gene location.** Two chains where only `t2`'s extra exon is antisense to a
  reference exon → the gene is `antisense_exonic`.
- **Hash collision.** Monkeypatched short hash → exit 1 naming both genes.
- **gffread.** It accepts `lnc_RNA` / `mRNA` transcripts with exon children. The existing
  validation already runs; assert that it still passes and that `gffread -E` sees one transcript
  per chain.

## Verification

1. **Unit tests.** `python3 tests/curate_models/test_curate_models.py`.
2. **Regression on `results_OctVul`** (all genes have `n_chains=1`). Run `bin/curate_models.py`
   directly on `targets/Octopus_vulgaris_6645/annotation/Octopus_vulgaris_6645.allModels.lnc_RNA.raw.gff3`
   with the samplesheet reference. Expect:
   - 54 genes, 54 transcripts;
   - `ref_location` tally unchanged (42 intergenic / 6 antisense_exonic / 3 antisense_intronic
     / 3 intronic);
   - the same representative model per gene;
   - every non-`novel` row of `curation.tsv` identical.
3. **Multi-chain behaviour on apiaceae** (`results_new`, e.g. `Ammi_trifoliatum_931225`: 3
   two-chain genes, from `annotation/*.allModels.lnc_RNA.raw.gff3`). Expect:
   - transcripts = genes + 3;
   - each 2-chain gene has `t1`/`t2` with distinct chains that share ≥ 1 SJ;
   - `t1` has ≥ species than `t2`;
   - the gene span covers both transcripts.
4. **Downstream readers.** gffread, gffcompare (`-r` the curated file vs. itself: 100 % Sn/Pr)
   and a quick IGV/JBrowse load accept `lnc_RNA` rows and the merged file.
5. **Pipeline.** `-profile test` (lncRNA only) and `-profile test --include_mrna --include_decoy`:
   CURATE_MODELS outputs publish under `targets/<t>/curated/`; the MultiQC curation-novel table
   shows `transcripts`; the run summary is unchanged.

## Not changed

- **Curation method.** Exact chains, ≥ `curate_min_species`, reference filter, SJ linking.
- **Upstream files.** The per-source projected GFF3s and `allModels` files keep their full
  alignment attributes; they are where the dropped metrics live.
- **Run-summary tables.**

## Resolved questions (2026-10-08)

- **mRNA track without CDS.** Accepted: curated `mRNA` rows have exons only. The track is a
  positive control, not part of the deliverable.
- **Source ID prefix.** Accepted: `representative_transcript` / `supporting_transcripts` carry
  the post-`RENAME_FASTA_HEADERS` tid, so Ensembl's `transcript:` prefix is stripped. The tid
  matches the projected GFF3s exactly, and the source GFF3 after re-adding the prefix.
- **Stale files from older runs.** Left in place; not the pipeline's concern.

## Results (2026-10-08)

- **Unit tests.** `tests/curate_models/test_curate_models.py` passes all 96 checks, including
  the new cases:
  - 20: longest representative;
  - 21: equal-support ordering by length;
  - 22: footprint location;
  - ID format, exon numbering and uniqueness;
  - gffread transcript count;
  - track-distinct hashes;
  - `--id-prefix`;
  - stability under reversed input and extra sources;
  - hash collision.
- **`results_OctVul` regression** (run directly against the samplesheet reference):
  - 54 genes and 54 transcripts;
  - every old representative is a `t1`, with identical exons and identical `ref_*` tags;
  - every non-`novel` row of `curation.tsv` is identical. The `novel` rows differ only by the
    new `transcripts` / `transcripts_per_gene_mean`.
- **apiaceae `Ammi_trifoliatum_931225`** (no reference):
  - 111 genes, 114 transcripts (3 two-chain genes);
  - every old representative is a `t1` with identical exons;
  - `genes_<n>sp` shifts (81/23/7 → 79/22/10) because it now counts species over all of a gene's
    chains;
  - transcript-level medians are now over 114 transcripts.
- **Downstream readers.**
  - gffread reads one transcript per `lnc_RNA` row (114 and 54);
  - gffcompare of the curated file against itself gives 100 % at the intron-chain and
    transcript levels.
