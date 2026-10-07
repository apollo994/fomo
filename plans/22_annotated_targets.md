# 22 — Annotated targets: a `target` may carry a gff3, used as the curation reference

## Status: implemented — merged to `main` @ `3dc0ce4` (squashed; detailed commit in tag `archive/plans-21-23`) — verification in **Results**

## Decisions (2026-10-07)

| Question | Decision |
|---|---|
| What a pure `target`'s gff3 is for | **Curation only** (plans/21). It is the reference that **(1)** removes curated models overlapping reference exons, **(2)** tags each curated gene with its nearest reference gene and its `ref_location`, and **(3)** is merged with the curated genes into `<target>.curated.<gtype>.merged.gff3.gz`. The target's raw GFF stats also appear in its report and in the run summary ("before" counts). |
| Benchmarking / top-N | **Unchanged.** Only `both` rows are benchmarked: the gffcompare reference is the species' own filtered ★ annotation (plans/20), which exists only for donors. A target with a gff3 is **not** benchmarked and gets no top-N consensus. To benchmark it, make it `both` (it then also donates). |
| `source` / `both` rows | Unchanged: `gff3` required. |
| Nearest-gene tag | On every curated gene: `ref_gene_id`, `ref_gene_name`, `ref_gene_biotype`, `ref_gene_orientation` (`sense` / `antisense` relative to the curated gene) and `ref_distance` (0 = overlapping, else the coordinate gap, adjacent = 1). It replaces plans/21's `ref_nearest` (type only). `ref_location` is unchanged. **Refined during implementation:** the tagged gene is the one the `ref_location` class refers to: the gene of the overlapping opposite-strand exon (`antisense_exonic`), of the innermost containing same-strand transcript (`intronic`), of an overlapping transcript (`sense_span_overlap` / `antisense_intronic`), or the head-to-head partner (`divergent`). Otherwise (`intergenic`) it is the nearest gene on either strand. Without this, an `antisense_exonic` gene could have been tagged with an overlapping *sense* gene, which contradicts its own class. |

## Context

- **Why the rule exists.** Since plans/20 a `target` row with a `gff3` aborts at launch
  (`assets/schema_input.json`, second `allOf` rule). The reason was benchmarking: its reference is
  the species' own filtered annotation, which only donors have. So "annotated target" meant
  `both`, and the only way to give a target a reference was to make it donate to every other
  target too.
- **Why that is now a problem.** Curation (plans/21) needs a reference for any target, and
  donating is not always wanted: a low-quality or out-of-clade annotation would pollute
  everyone's `allModels`, and the self-projection would add to the target's own report.
- **The pipeline is already almost ready:**
  - `CURATION` reads the reference from the `gff3` slot of `ch_targets`, whatever the role
    (`subworkflows/local/curation.nf`).
  - `BENCHMARKING` takes its references from `PREPROCESSING.out.filtered_gff3` (donors only). A
    target-with-gff3 has none, so the projection ⋈ reference inner join drops it. There is no
    `failOnMismatch` on that path.
  - `BENCHMARKING` already runs `GFF_STATS_TARGET` (`kind: 'raw'`) for every `ch_targets` row
    with a gff3, so the target's annotation stats appear in its report without any change.
  - `PREPROCESSING`, decoys and `CONSENSUS_TOP` only see `ch_sources` or benchmarking output.
    Nothing new reaches them.

## Changes

### 1. Samplesheet schema (`assets/schema_input.json`)
- **Remove the second `allOf` rule** (`target` ⇒ `not required gff3`). Keep the first
  (`source` / `both` ⇒ `gff3` required).
- **Reword the `role` `errorMessage`** and the first rule's message:
  - `target` = projected onto; it does not donate. Its optional `gff3` is used only as the
    curation reference and is not benchmarked.
  - To benchmark against its own annotation, use `both`, which also makes it a donor.

### 2. Workflow (`workflows/fomo.nf`)
- **Launch log line.** When there is at least one target with a gff3, log one info line listing
  the annotated non-donor targets: "used as curation reference only; not benchmarked (use role
  'both' to benchmark)". This makes the non-benchmarking visible, so it is not a surprise.
- **No change otherwise.** `ch_roles` already writes `has_gff3` per species.

### 3. Run summary (`bin/run_summary_tables.py`)
- **"Targets benchmarked (with GFF3)"** currently counts `is_target and has_gff3`. It becomes
  `role == 'both'`, which is exactly what is benchmarked.
- **New row** "Targets with a reference annotation" (`is_target and has_gff3`), and a matching
  `n_targets_with_reference` in `run_summary.json`. `n_targets_annotated` keeps its current
  meaning (benchmarked targets), so existing consumers of the JSON see no change.
- **"… role target (annotated only)"** becomes "… role target (projected onto, not a donor)".
- **Per-target "before" counts already use `has_gff3`.** An annotated pure target now shows its
  existing lncRNA count, which is correct: it is a real annotation. Check that the ranking and
  top-N sections do not assume `has_gff3 ⇒ benchmarked` anywhere else (grep `has_gff3`,
  `annotated`).

### 4. Nearest reference gene (`bin/curate_models.py`)
- **Gene identity per reference transcript span**, resolved while reading the reference:
  - **GFF3:**
    - `gene_id` = the transcript's `Parent` when that parent is a gene-level feature (`gene`,
      `pseudogene`, `ncRNA_gene`, or any feature that is a Parent and has no Parent itself).
    - An exon parented directly by a gene → that gene.
    - A transcript without a parent → the transcript itself.
    - `gene_name` = the gene's `Name` → `gene_name` → `gene` attribute, then the transcript's
      `gene` (RefSeq), else `.`.
    - `gene_biotype` = the gene's `gene_biotype` → `biotype` → `gene_type`, else the
      transcript's type label (e.g. `lnc_RNA`, `transcript:misc_RNA`).
  - **GTF:** `gene_id`, `gene_name`, `gene_biotype` / `gene_type` from the exon or transcript
    attributes.
  - Stored as one entry per gene (`gene_id → (name, biotype, strand)`). The span index keeps a
    small integer label per transcript span instead of the type string. Memory must stay within
    plans/21's 772 MB worst case (Mercenaria); re-measure.
- **Nearest gene.** Minimum distance over the transcript spans of both strands (the existing
  `RefIndex.nearest`, extended to return the label). Ties: overlapping before distant, then
  sense before antisense, then smaller gene ID. Distance is 0 for any span overlap, including
  intronic and antisense. For `divergent` / `intergenic` it is the gap in bp.
- **Attributes on the representative transcript** (and `ref_gene_id` + `ref_distance` on the
  gene row too):
  - `ref_gene_id`, `ref_gene_name`, `ref_gene_biotype`, `ref_gene_orientation`,
    `ref_distance`;
  - all `.` when the seqid has no reference transcript; absent without `--ref`;
  - `ref_nearest` is removed.
- **Report** (`novel` table): add `intergenic_distance_{min,p25,median,mean,p75,max}` over
  `ref_location=intergenic` genes. `bin/curation_to_mqc.py` adds `median_intergenic_distance`
  to the novel MultiQC table.
- **Consistency guard (exit 1 if violated):**
  - `ref_distance == 0` for every overlapping class (`antisense_exonic`, `intronic`,
    `sense_span_overlap`, `antisense_intronic`), and `> 0` for `divergent` / `intergenic`;
  - `ref_gene_orientation == antisense` for `antisense_*`;
  - `ref_gene_orientation == sense` for `intronic`.

### 5. Tests
- **Fixtures** (`tests/curate_models/test_curate_models.py`):
  1. gene `Name` vs RefSeq `gene=` vs no name (`.`);
  2. `gene_biotype` vs fallback to transcript type;
  3. GTF `gene_name` / `gene_type`;
  4. exon parented directly by a gene;
  5. nearest on the opposite strand closer than on the same strand → `antisense`;
  6. distance tie → sense wins;
  7. seqid without reference → `.`;
  8. the consistency guard per `ref_location` class;
  9. intergenic distance stats.
- **Default test samplesheet** (`assets/samplesheet.csv`):
  - Give the pure target Lycaena hippothoe its gff3 (the file is already in `assets/test_data/`).
    The default `-profile test` then exercises the new path: 2 sources (alciphron, thersamon)
    plus a real reference, so the reference filter actually removes models.
  - Keep the reference-less path as `assets/samplesheet_target_noref.csv` (today's sheet).
- **Schema checks:** a `target` with a gff3 validates; a `source` without one still aborts with
  the reworded message.

## Files to touch

| File | Change |
|---|---|
| `assets/schema_input.json` | drop the target-without-gff3 rule; reword messages |
| `assets/samplesheet.csv`, `assets/samplesheet_target_noref.csv` (new) | hippothoe gets its gff3; keep the old sheet |
| `workflows/fomo.nf` | launch info line for annotated non-donor targets |
| `bin/run_summary_tables.py` | benchmarked = `both`; new "with a reference annotation" row / JSON key; role label |
| `bin/curate_models.py` | gene identity in the reference index; nearest gene tag; intergenic distance; guard |
| `bin/curation_to_mqc.py` | `median_intergenic_distance` column |
| `assets/multiqc/sections.yml` | novel-table description mentions the new columns |
| `tests/curate_models/test_curate_models.py` | fixtures above |
| `CLAUDE.md`, `README.md` | Samplesheet & roles section; "annotated targets are exactly the `both` rows" statements (CLAUDE.md, benchmarking.nf comments) |
| `plans/20_longest_isoform.md` | one-line note under its Decisions: the "pure targets with a gff3 not allowed" decision is superseded by plans/22 for curation |
| `subworkflows/local/benchmarking.nf` | comment only: "a gff3 on a pure target is rejected" → "is used for curation, not benchmarking; dropped by the inner join" |

## Verification

1. Fixtures pass on the host and in the `CURATE_MODELS` container (gffcompare 0.12.6).
2. **Real data, CLI**, on the 10 plans/21 rev-2 targets:
   - curated gene sets and `ref_location` are identical to the current outputs (only the new
     attributes are added);
   - the consistency guard holds;
   - report the distribution of nearest-gene biotypes and intergenic distances;
   - peak RSS on Mercenaria is ≤ plans/21's 772 MB, or the change is explained.
3. **`-profile test,crg`** with the new default sheet:
   - hippothoe has a merged file, raw GFF stats in its report and curation genes tagged with
     reference genes;
   - its `discarded_by_reference` table is non-zero;
   - no gffcompare / top3 / select_top_sources outputs for hippothoe;
   - the run summary counts it under "with a reference annotation" and not under
     "benchmarked";
   - no `targets/null/`.
4. **`--input assets/samplesheet_target_noref.csv`:** identical to the plans/21 phase B test
   (reference-less path).
5. **Schema:** `target` + gff3 accepted; `source` without gff3 rejected with the new message.
6. `nextflow lint .`: no new errors.

## Out of scope

- **Benchmarking annotated non-donor targets.** It would need FILTER_TRANSCRIPT / TD2 /
  FILTER_ANNOTATION run on targets that do not donate, and a change to the plans/20
  filtered-annotation contract. Possible follow-up if wanted.
- **A separate curation reference column** (e.g. a fuller annotation than the one donated by a
  `both` row). `gff3` serves both roles for `both`.

## Results (2026-10-07)

**Fixtures:** 75 checks pass, on the host and in the `CURATE_MODELS` container (gffcompare 0.12.6).
- New case 19 covers:
  - the gene `Name` vs RefSeq `gene=` vs no name;
  - `gene_biotype` vs the fallback to the transcript type;
  - GTF `gene_name` / `gene_type`;
  - exons parented directly by a gene;
  - the opposite-strand gene being nearer → `antisense`;
  - a distance tie → sense wins;
  - a seqid without reference → `.`;
  - the gene row carrying `ref_gene_id` + `ref_distance`;
  - intergenic distance stats;
  - the consistency guard rejecting three contradictory tags.
- Case 18 now checks the divergent partner (`rna-dv`) and the nearest gene (`rna-far`).

**Real data, CLI, HEAD vs new.** The 10 plans/21 rev-2 targets plus Mercenaria:
- curated gene sets, `ref_location` and every non-tag attribute are **identical** on all 11;
- `ref_distance` differs for 10 genes, **all `divergent`**: they are now tagged with their
  head-to-head partner rather than the nearest gene (e.g. 219 → 577 bp);
- across 5,470 genes, the tagged gene's biotype is:

  | biotype | genes |
  |---|---|
  | protein_coding | 3,604 |
  | lncRNA | 966 |
  | misc_RNA | 391 |
  | Y_RNA | 350 |
  | pseudogene | 49 |
  | tRNA | 45 |
  | `.` (seqid without reference) | 30 |

- 3,082 of them have a gene name;
- orientation is consistent with every class;
- intergenic distance: median 5.3 kb (p25 1.3 kb, p75 20.7 kb);
- Mercenaria peak RSS 791 → 931 MB, from the per-gene name/biotype tables. Still well inside
  `process_low` (4 GB), so accepted.

**Schema:**
- a `source` without a gff3 is rejected with the reworded message;
- a `target` with a gff3 validates, and the launch log lists it as "curation only (not
  benchmarked)".

**`-profile test,crg`, new default sheet** (hippothoe `target` + gff3), not resumed: 79 tasks
succeeded.
- Hippothoe:
  - 17 supported chains; 14 lost to the reference (primary type lnc_RNA 9, mRNA 5);
  - 3 curated genes, one `antisense_exonic` and two `intergenic`, tagged e.g. METTL3 (404 bp) and
    COBLL1 (4.5 kb);
  - a merged file;
  - raw target GFF stats in its report;
  - **no** gffcompare / select_top_sources / top3 outputs.
- Run summary: "Targets with a reference annotation (curation)" 2, "Targets benchmarked (role
  both)" 1; hippothoe's species row has `has_gff3=yes`, `benchmarked=no`. No `targets/null/`.

**`--input assets/samplesheet_target_noref.csv`** (resumed): hippothoe is curated
reference-free (17 genes, `no_reference`, no merged file).

**`nextflow lint .`:** the same 7 errors as before, none new.
