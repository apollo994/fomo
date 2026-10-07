# 20 — One filtered annotation per species: longest isoform, used for projection, reference and results

## Status: implemented (2026-10-01), verified against baseline on `-profile test`

Implemented on branch `longest-isoform` as planned; one detail changed: the funnel row's
`Sample` is `<species>.<feature_type>`, not the bare species — with `--include_mrna` a
species has two rows and MultiQC keys table rows by sample name (`run_summary_tables.py`
strips the known `.<feature>` suffix by length).

**Verification performed** (real `-profile test,singularity` runs, 1-CPU cap, baseline =
`main` @ `5560ac9` in a separate worktree):

- `filter_transcript.sh`: without `--longest`, byte-identical to `main` on all 5 test GFF3s
  × {lnc_RNA, mRNA}; mawk (`ubuntu:22.04`) and gawk outputs identical; picks match an
  independent Python re-derivation on all 10 (incl. every tie: 0–1 lncRNA, 3–11 mRNA per
  species); reruns byte-identical; fixture covers tie, mono-exonic longest, orphans,
  `Parent=a,b`, mixed mRNA/lnc_RNA gene, `transcript:foo:bar`, empty input.
- `filter_gff_by_id.py`: keep mode == drop mode on real data (phlaeas, every 5th dropped);
  all three exit-1 checks fire; drop mode byte-identical to `main` on the 4 real
  `FILTER_ALLMODELS` / `FILTER_LNC_GFF` inputs of the baseline run.
- `nextflow lint .`: the same 7 errors as `main`, all in untouched files.
- `--longest_isoform false` vs baseline: `spliced_fasta` and `<sp>.lnc_RNA.filtered.gff3`
  byte-identical for both sources; SeqKit and TD2 tables identical. Only the target-side
  `filtered` stats row differs (alciphron reference 67 → 64 lncRNA: now TD2-filtered), and
  with it transcript Sn (self 88.1 → 92.2, thersamon 26.9 → 28.1); Pr unchanged.
- Default run vs baseline (target alciphron, lncRNA, transcript level Sn | Pr):
  self 88.1 | 98.3 → 91.5 | 98.2; thersamon 26.9 | 28.1 → 28.8 | 28.3; missed loci
  7/62 → 4/59. Top-source ranking unchanged (one non-self donor in the test sheet).
- Funnel: alciphron 67 = 0 single-exon + 5 non-longest + 3 coding + 59 kept; thersamon
  85 = 0 + 8 + 3 + 74; no consistency warning; `transcript_filter` renders once per report.
- Chain of custody: the reference staged for gffcompare has the same md5 as
  `sources/Lycaena_alciphron_282377/annotation/…lnc_RNA.filtered.gff3`; projected FASTA
  records == filtered GFF3 transcripts (via the id map) for both sources.
- `--include_mrna --include_decoy`: mRNA filtered GFF3 published, 1 transcript/gene (887,
  864); decoys built from the pre-TD2 candidates (lncRNA 62 / 77 = n_selected).
- No `targets/null/` or `sources/null/`; `sources/` holds exactly the two `source`/`both`
  rows. Samplesheet: a `target` row with a `gff3` aborts at launch with the new message;
  the committed sheet validates. `assets/samplesheet_full.csv` (phlaeas → `both`) was **not**
  validated: its `assets/annotation_downloads/` files are not present on this machine.

## Decisions (2026-10-01)

| Question | Decision |
|---|---|
| Representative transcript | longest spliced isoform per gene, by summed exon length; AGAT-free, inside `filter_transcript.sh` |
| Tie-break | lexicographically smallest transcript ID (generic, not `Ensembl_canonical`) |
| One file | ★ `<sp>.<ft>.filtered.gff3` is built **from** the projected FASTA via an id map written at rename time — one extraction, no re-extraction |
| gffcompare reference | ★ — `FILTER_TARGET` deleted; the reference is longest-isoform **and** TD2-filtered |
| Pure targets with a `gff3` | **not allowed** — a `target` row with a `gff3` fails samplesheet validation (§4). Only `source`/`both` rows carry an annotation, so ★ exists exactly for the sources and every annotated target is a `both`. *Superseded by plans/22 (2026-10-07): a pure target may carry a `gff3`, used as its curation reference only; benchmarking stays `both`-only.* |
| Decoys | relocated from `FILTER_TRANSCRIPT` output (spliced + longest, pre-TD2) — TD2 is irrelevant to a decoy |
| Results | `sources/<sp>/annotation/` — one per `source`/`both` row |
| Funnel reporting | new `dropped_non_longest` column, fed by a `FILTER_TRANSCRIPT` counts table |

## Context

Two things are wrong with how a species' annotation flows through the pipeline today.

**Every isoform is transferred.** `FILTER_TRANSCRIPT` keeps every multi-exon transcript of
the enabled type, so a gene with five isoforms contributes five near-redundant queries. We
want **one representative transcript per gene — the longest**. AGAT
(`agat_sp_keep_longest_isoform.pl`) is not an option: the pipeline dropped AGAT entirely
(plans/10, plans/11) and must not re-acquire it.

**There is no single "filtered annotation".** The same samplesheet GFF3 is filtered by
separate paths that stop at different points, and are kept in sync by luck:

| Consumer | Built from | TD2-filtered? |
|---|---|---|
| projected FASTA (`spliced_fasta`) | `FILTER_TRANSCRIPT` → gffread → rename → TD2 kept FASTA | yes |
| `filtered_gff3` / input "filtered" stats | `FILTER_TRANSCRIPT` → `FILTER_LNC_GFF` (TD2 coding ids) | yes |
| gffcompare reference | `BENCHMARKING:FILTER_TARGET`, a second, independent run of the same filter | **no** |
| decoys | `FILTER_TRANSCRIPT` (pre-TD2) | no — and that is correct, see §6 |

The FASTA and the filtered GFF3 are tied together only by **two independent id
normalisations that do not agree**:

- `RENAME_FASTA_HEADERS` turns a header into a bare tid by `split($1, a, ":")` and taking
  `a[2]` — i.e. *everything between the first and second colon*;
- `TD2_CODING_FILTER --stage input` cuts the renamed header back at the first `|`;
- `filter_gff_by_id.py` strips a fixed prefix list (`transcript:`, `gene:`, `rna-`, …) from
  the GFF side and the drop list before comparing.

For Ensembl (`transcript:ENSGXNT…`) and NCBI (`rna-XM_…`, no colon) these happen to meet.
For an id with a second colon (`transcript:foo:bar` → rename gives `foo`, the GFF side
gives `foo:bar`) or an unlisted prefix, a coding transcript is dropped from the FASTA but
**silently kept** in the GFF3 — and nothing checks.

**Target design:** one file per source species (`source` | `both`) and feature type,
`<species>.<ft>.filtered.gff3` (★) — spliced, longest isoform per gene, lncRNA
TD2-noncoding — that

1. **describes exactly the sequences sent to projection**, record for record, checked at
   runtime;
2. **is the gffcompare reference** for projections onto that species;
3. **is published** in the results.

Profile of the five test GFF3s (Ensembl-style `gene:` → `transcript:` → `exon`), counted on
the raw annotation:

| track | genes | multi-isoform genes | spliced-length ties | longest ≠ `Ensembl_canonical` |
|---|---|---|---|---|
| lnc_RNA | 62–112 | 4–15 | 0–1 | 0 |
| mRNA | 889–930 | 326–387 | 4–12 | 70–106 (~8–11 %) |

So `-profile test` exercises the change (and the tie-break) on lncRNA alone, and heavily
with `--include_mrna`.

## New flow

```
samplesheet gff3 — every SOURCE (role source | both; the only rows allowed a gff3, §4)
 └─ FILTER_TRANSCRIPT --longest           <sp>.<ft>.gff3   (candidates: spliced + longest)
     │
     ├─ EXTRACT_SEQUENCES (gffread -w)    ONE extraction, as today
     │   └─ RENAME_FASTA_HEADERS          <sp>.<ft>.spliced.renamed.fasta
     │                                    <sp>.<ft>.id_map.tsv   (renamed header → GFF ID)  NEW
     │       ├─ lnc_RNA: TD2_NONCODING    → kept (non-coding) FASTA
     │       └─ mRNA:    as is
     │                │
     │                ├──────────────────────────────▶ spliced_fasta → PROJECTION
     │                │                                 does NOT wait for ★
     │                ▼
     ├────────▶ FILTER_ANNOTATION(candidates gff3, id_map, final FASTA)
     │            keep exactly the transcripts whose header is in the FASTA
     │            ★ <sp>.<ft>.filtered.gff3
     │               ├─ BENCHMARKING reference (role both) → GFFCOMPARE_BATCH,
     │               │    target GFF stats, target_refs → CONSENSUS_TOP
     │               ├─ input "filtered" GFF stats
     │               └─ published: sources/<sp>/annotation/
     │
     └─ RELOCATE_LOCI → EXTRACT_DECOY_SEQUENCES   (--include_decoy; unchanged input)
```

One extraction, one rename, one TD2 pass per species and feature type — the same work as
today. ★ is derived **from** the projected FASTA, so the two cannot disagree; projection
runs in parallel with building ★.

## Design

### 1. Longest-isoform selection inside `bin/filter_transcript.sh` (no new process)

The script is a two-pass awk: pass 1 collects candidates of `TTYPE` and counts exons per
parent; a boundary block promotes multi-exon candidates into `tx[]` / `parent_of_tx[]`;
pass 2 prints the lineage of everything in `tx[]`. Selection slots into that structure:

- **Pass 1:** alongside `exon_count[p[j]]++`, accumulate `exon_len[p[j]] += $5 - $4 + 1`.
  The loop already iterates every id of a comma-separated `Parent=`, so an exon shared by
  several transcripts credits each.
- **Boundary block:** for each candidate with `exon_count >= 2`, group key = its full
  `Parent=` string (its own id when it has none — an orphan competes only with itself).
  Per group keep the candidate with the greatest `exon_len`; **ties go to the
  lexicographically smallest transcript ID** (awk `for (k in a)` order is unspecified, so
  the rule must be explicit for reproducibility). Build `tx[]` / `parent_of_tx[]` from the
  **winners only**. Pass 2 is untouched.
- **Flag:** `filter_transcript.sh [--longest] input.gff TRANSCRIPT_TYPE`, parsed with the
  `while/case` loop of `bin/bam_to_gff.sh` (unknown `-*` → exit 1). Without `--longest` the
  boundary block behaves exactly as today.

| Decision | Why |
|---|---|
| Length = **sum of exon lengths**, not span, not CDS | It is exactly what `gffread -w` extracts and minimap2 aligns; span is intron-inflated. AGAT ranks mRNA by CDS, but we project cDNA, so exon length is consistent across both tracks. |
| Select **after** the multi-exon filter | Otherwise a gene whose longest isoform is mono-exonic would lose the gene despite having a spliced isoform. |
| Select **before** TD2 | A gene whose longest lncRNA isoform is coding is dropped — a defensible call on the locus — and TD2 gets a smaller FASTA. Selecting after TD2 would need TD2 on every isoform plus a FASTA→GFF round-trip per gene. |
| Per (gene, feature type) | The script runs once per type, so a gene with both an `mRNA` and a `lnc_RNA` child yields one winner per track. |
| No new process | A separate `KEEP_LONGEST_ISOFORM` would add `F·S` short SLURM jobs — the overhead plans/17 removed. |

Not `gffread`: it has no keep-longest mode; emulating one (`--table @id,@geneid,@covlen` →
keep-list → filter) is two steps for the same result. Move to a Python script only if the
rule grows (e.g. CDS-length ranking for mRNA).

One param, `--longest_isoform` (default `true`), on the one place it applies:

```groovy
withName: 'FOMO:PREPROCESSING:FILTER_TRANSCRIPT' {
    ext.prefix = { "${meta.id}.${meta.feature_type}" }
    ext.args   = params.longest_isoform ? '--longest' : ''   // precedent: BAM_TO_GFF
}
```

### 2. `RENAME_FASTA_HEADERS` writes the id map

The rename step is the only place that sees both names, so it records them. Alongside the
renamed FASTA it emits `<prefix>.id_map.tsv`, one row per record, no header:

```
<renamed header token>\t<original gffread seqname>
ENSGXNT00000001570|lncRNA|Lycaena_alciphron_282377	transcript:ENSGXNT00000001570
```

`gffread -w` names each record by its transcript `ID=` verbatim, so column 2 is the exact
GFF3 transcript id — no normalisation anywhere. In the existing awk this is one extra
`print ... > map` next to the header rewrite; the module gains
`tuple val(meta), path("*.id_map.tsv"), emit: id_map`. The renamed FASTA is byte-identical
to today. Decoy FASTAs get a map too; it is simply not consumed.

The `a[2]` colon split itself stays as is: it only shapes the *renamed* header, whose
uniqueness `FILTER_ANNOTATION` now checks (§3), and changing it would change every
projected id downstream.

### 3. `FILTER_ANNOTATION` builds ★ from the FASTA that is actually projected

New local module `modules/local/filter_annotation.nf`, `process_single`, `python:3.11`:

```groovy
input:
tuple val(meta), path(gff3, stageAs: 'input/*'), path(id_map), path(fasta)
output:
tuple val(meta), path("*.filtered.gff3"), emit: gff3
script:
filter_gff_by_id.py --gff ${gff3} --keep-fasta ${fasta} --id-map ${id_map} --output ${prefix}.gff3
```

`bin/filter_gff_by_id.py` gains a **keep mode** (`--keep-fasta` + `--id-map`, mutually
exclusive with `--drop-ids`):

1. read the map; a duplicate renamed header → **exit 1** naming it;
2. read the FASTA record names; one not in the map → **exit 1**;
3. translate to GFF ids and keep exactly those transcripts (`drop = all transcript ids −
   keep`, then the existing lineage logic: children follow their transcript, a gene
   survives iff ≥ 1 transcript child survives). Matching is **exact** — `normalise()` is
   not used in keep mode;
4. a kept id not present in the GFF3 as a transcript → **exit 1**;
5. so on success, transcripts in ★ == records in the projected FASTA, by construction.

Drop mode (`--drop-ids`, used by `PROJECTION:FILTER_ALLMODELS`) is unchanged.

Wiring in `preprocessing.nf` — the "final FASTA" is the one already computed today as
`ch_renamed_final` (TD2 kept for lncRNA, renamed as-is for mRNA), real track only:

```groovy
ch_final_fasta = ch_renamed_final.filter { m, _fa -> !m.decoy }

// [meta.id, feature_type] is the key: one candidate GFF3, one id map and one final
// FASTA per species and feature type.
// failOnMismatch: a species whose FASTA or map never arrives must abort the run, not
// silently vanish from the results, the reference and the report.
ch_annot_in = FILTER_TRANSCRIPT.out.gff3
    .map { m, g  -> tuple(m.id, m.feature_type, m, g) }
    .join(RENAME_FASTA_HEADERS.out.id_map.filter { m, _t -> !m.decoy }
              .map { m, t -> tuple(m.id, m.feature_type, t) }, by: [0, 1], failOnMismatch: true)
    .join(ch_final_fasta.map { m, fa -> tuple(m.id, m.feature_type, fa) },
          by: [0, 1], failOnMismatch: true)
    .map { _id, _ft, m, g, map, fa -> tuple(m, g, map, fa) }

FILTER_ANNOTATION(ch_annot_in)
ch_filtered_gff3 = FILTER_ANNOTATION.out.gff3          // ★
```

This **replaces** `FILTER_LNC_GFF` (and `TD2_NONCODING.out.coding_ids` loses its only
consumer — keep emitting it, it is cheap and useful when debugging) and needs no
`empty_ids.txt` trick for mRNA: mRNA's final FASTA is simply every extracted record.
Output name `<sp>.<ft>.filtered.gff3` — the lncRNA filename is unchanged.
`filter_gff_by_id.py` preserves line order, so ★ stays grouped by seqid
(gff-feature-stats requirement).

Side benefit: if gffread ever skips a transcript (e.g. a seqid missing from the FASTA), it
is now absent from ★ too, instead of being in the reference/stats but never projected.

### 4. Samplesheet: a pure `target` row may not carry a `gff3`

For ★ to be the only reference, every annotated target must also have been through
PREPROCESSING — which today runs only for sources. Rather than widening PREPROCESSING to
pure targets (a role split inside it, source-only filters on `spliced_fasta`, the decoys
and the input-side MultiQC rows, and a real risk of a pure target leaking into PROJECTION
as a donor), the samplesheet forbids the case (decided 2026-10-01): **a `gff3` belongs to a
species that donates it.**

`assets/schema_input.json` — a second `allOf` entry next to the existing
"source/both MUST have a gff3" rule:

```json
{
    "if": {
        "properties": { "role": { "const": "target" } },
        "required": ["role"]
    },
    "then": {
        "not": { "required": ["gff3"] },
        "errorMessage": "A row with role 'target' must leave gff3 empty: a target is only projected onto. To also benchmark it against its own annotation, use role 'both' — note that this makes it a donor to every other target too."
    }
}
```

An **error, not a silent ignore**: a user who supplies a target `gff3` expects
gffcompare numbers, and dropping the file quietly would give them a run without
benchmarking and no hint why. The existing rule's `errorMessage` ("Only a pure 'target'
may leave it empty — it is then annotated but not benchmarked") stays correct; the
`role` enum's message ("'target' (annotated only)") gains "(no gff3)".

This relies on nf-schema treating an **empty CSV cell as an absent key** — which the
current test sheet already depends on (`Lycaena_hippothoe_580924` has an empty `gff3` and
passes the `gff3` `pattern`). Verify the new rule both ways (Verification §4).

What this buys:

- **PREPROCESSING is unchanged in scope**: it keeps taking `ch_sources`, `S` tasks, no
  role logic inside, nothing new that could leak into `spliced_fasta`.
- `workflows/fomo.nf` changes only the BENCHMARKING call:
  `BENCHMARKING(ch_targets, PROJECTION.out.gff3, PREPROCESSING.out.filtered_gff3)`.
- `T_g` (targets with a `gff3`) is now exactly the `both` rows.

Accepted consequence: **a species can no longer be benchmarked without also donating.**
Its own annotation then joins the merged alignment onto every other target and competes in
their top-N rankings. If "benchmark-only" targets are needed later, they come back as a
separate feature (e.g. a `reference` role that runs PREPROCESSING but is filtered out of
`spliced_fasta`) — not by re-allowing a `gff3` on `target`.

Also update `assets/samplesheet_full.csv`: `Lycaena_phlaeas_282391` is a `target` **with**
a `gff3` and would now fail validation — either change its role to `both` or blank its
`gff3` (owner's call).

### 5. BENCHMARKING takes its reference from ★ — `FILTER_TARGET` is deleted

```groovy
workflow BENCHMARKING {
    take:
    ch_target          // [ meta, fasta, gff3 ] — raw gff3 still used for kind:'raw' stats
    ch_projected_gff3
    ch_annotation      // [ meta(id, role, feature_type, decoy:false), filtered.gff3 ] × F·S
    main:
    // "is this species a target". Only 'both' can match (a pure target has no gff3,
    // §4), but the condition states the intent and survives a future 'reference' role.
    ch_target_refs = ch_annotation.filter { m, _g -> m.role in ['target', 'both'] }
    ...
```

- Drop `FILTER_TARGET` and its `modules.config` block; every `FILTER_TARGET.out.gff3` use
  (the paired join, the `kind: 'filtered'` target stats, the `target_refs` emit) reads
  `ch_target_refs`. Meta is shape-compatible: both joins (`benchmarking.nf`,
  `consensus_top.nf:104`) key on `[meta.id, meta.feature_type]`.
- Drop the `feature_types` take: references now arrive already split by the very list
  PREPROCESSING used, so the "two sides must agree" invariant in `fomo.nf` holds by
  construction. Update that comment.
- Pure targets drop out exactly as today's un-annotated ones: they are not sources, so
  they have no ★, and the inner join removes their projections.

Consequences — accepted, to be written into CLAUDE.md:

- **The reference is now TD2-filtered.** This reverses the current `FILTER_TARGET`
  comment ("this reference copy must not be"). The query side is TD2-filtered too
  (PROJECTION's projected-lncRNA pass), so scoring it against an unfiltered reference
  counted coding lncRNA the pipeline would never report as false negatives. Both sides are
  now "non-coding lncRNA".
- **Transcript-level Sn ≈ per gene.** The denominator drops from "every spliced isoform" to
  "one per gene".
- **Transcript-level precision can drop**: a source's longest isoform matching a
  *non-longest* isoform of the target's ortholog used to be a TP, now it is an FP (and the
  reference's longest a FN). Mostly an mRNA effect; multi-isoform lncRNA genes are rare.
- **Self-pairs become an exact ceiling.** For a `both` species the reference *is* the
  description of the FASTA it projected, so X→X loses only what alignment and the
  projected TD2 pass lose.
- **Numbers are not comparable to pre-plan-20 runs**; the top-N ranking may shift.

### 6. Decoys stay on `FILTER_TRANSCRIPT` output

`RELOCATE_LOCI` keeps reading `FILTER_TRANSCRIPT.out.gff3` — its current input — restricted
to sources. TD2 is irrelevant to a decoy: only the model's *shape* (exon count, exon and
intron lengths) is borrowed, and the sequence is extracted from intergenic DNA at the new
site, so the original locus' coding potential says nothing about it. Not waiting on ★ also
keeps the decoy branch off TD2's critical path.

Decoys do pick up the longest-isoform selection, because it happens inside
`FILTER_TRANSCRIPT`. That is a fix in its own right: today every isoform is relocated as an
independent model, so multi-isoform genes are over-represented in the decoy set. Expect
different decoy numbers with `--include_decoy` (the `decoy_cap` sample is drawn from a
smaller pool), but no wiring change. Update the existing `preprocessing.nf` comment so it
says *why* decoys skip TD2, instead of "remain an unfiltered null baseline".

### 7. Publish ★

Published under `sources/` (decided 2026-10-01):

```
sources/<species>/annotation/  <species>.lnc_RNA.filtered.gff3
                               <species>.mRNA.filtered.gff3      (--include_mrna)
```

One entry per `source`/`both` row — exactly the species that donate (§4). A sibling of
`targets/` and `summary/`; a `both` species appears under both `sources/` and `targets/`,
which hold different things (its input vs. what was projected onto it).

```groovy
// ── sources/ — each source's filtered annotation, a sibling of targets/ ──
// Keyed on meta.id: target-independent (PREPROCESSING runs once per species), so
// exactly one copy. Describes exactly what was projected, and is the gffcompare
// reference for projections onto this species.
withName: '.*:FILTER_ANNOTATION$' {
    publishDir = [
        path:   { "${params.outdir}/sources/${meta.id}/annotation" },
        mode:   params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Fix that section's header comment ("Every deliverable is target-scoped") and rule 1.

### 8. Report the isoform drop separately (run-summary source funnel)

`section_source_funnel` (`bin/run_summary_tables.py:395`) computes
`dropped_single_exon = raw_lncRNA_transcripts − TD2 n_in`. With longest-isoform selection
TD2's input is spliced **and** longest, so that column would silently become "single-exon
+ non-longest isoform" under the old label. Nothing upstream records the count between the
two filters, so `FILTER_TRANSCRIPT` reports its own funnel (decided 2026-10-01: new column,
not a relabel).

- **`bin/filter_transcript.sh --counts FILE --sample NAME`**: at the boundary block (where
  all three numbers are known) append one row to `FILE`:

  ```
  Sample	feature	n_transcripts	n_spliced	n_selected
  Lycaena_alciphron_282377	lnc_RNA	67	67	62
  ```

  `n_transcripts` = candidates of `TTYPE`, `n_spliced` = with ≥ 2 exons, `n_selected` =
  after longest-isoform selection (`== n_spliced` without `--longest`).
- **`modules/local/filter_transcript.nf`**: pass
  `--counts ${prefix}_transcript_filter_mqc.tsv --sample ${meta.id}`; new emit
  `tuple val(meta), path("*_transcript_filter_mqc.tsv"), emit: mqc`.
- **`preprocessing.nf`**: mix `FILTER_TRANSCRIPT.out.mqc` into `ch_mqc_files` (every
  PREPROCESSING row is a source, so no role filter is needed).
- **MultiQC** (per CLAUDE.md "Adding a new custom-content section"): unique suffix
  `*_transcript_filter_mqc.tsv`, a `transcript_filter` `custom_data` block + `sp` pattern in
  `assets/multiqc/sections.yml`, a `report_section_order` slot just before `td2_coding` in
  `main.yml`, and the suffix in `extra_fn_clean_exts`. No `# id:` header in the file
  (rule 1). It is a source-side table, broadcast into every per-target report like TD2's.
- **`run_summary_tables.py`**: `self.filter_rows = read_tsvs(pick("_transcript_filter_mqc.tsv"))`;
  the funnel becomes

  | column | from |
  |---|---|
  | `dropped_single_exon` | `n_transcripts − n_spliced` |
  | `dropped_non_longest` **(new)** | `n_spliced − n_selected` |
  | `dropped_coding` | TD2 `n_coding` (unchanged) |

  plus two consistency warnings in the style of the existing one at :418:
  `n_selected ≠ TD2 n_in` and `n_transcripts ≠` raw-GFF-stats lncRNA count. Add the new
  column to the funnel bargraph (`run_source_funnel_mqc.json`) and to `run_summary.json`.
  Update the docstring's `run_source_funnel` line ("raw → spliced → longest → non-coding").

## Task counts

| | today | after |
|---|---|---|
| `FILTER_TRANSCRIPT`, `EXTRACT_SEQUENCES`, `RENAME_FASTA_HEADERS`, `TD2_NONCODING` | unchanged | unchanged |
| `FILTER_LNC_GFF` → `FILTER_ANNOTATION` | `S` | `F·S` (`+S` only with `--include_mrna`) |
| `BENCHMARKING:FILTER_TARGET` | `F·T_g` | — |

Net, lncRNA only: **−T_g** tasks; everything else is unchanged.

## Files to touch

| File | Change |
|---|---|
| `bin/filter_transcript.sh` | `--longest`, `exon_len` tally, per-gene selection; `--counts`/`--sample` funnel row (§8); header documents rule, tie-break, order vs. multi-exon filter |
| `modules/local/filter_transcript.nf` | pass `${task.ext.args ?: ''}` and `--counts`/`--sample`; new `mqc` emit |
| `assets/multiqc/sections.yml`, `assets/multiqc/main.yml` | `transcript_filter` section: `custom_data` + `sp`, `report_section_order`, `extra_fn_clean_exts` (§8) |
| `modules/local/rename_fasta_headers.nf` | write + emit `*.id_map.tsv` |
| `bin/filter_gff_by_id.py` | keep mode: `--keep-fasta` + `--id-map`, exact matching, three exit-1 checks |
| `modules/local/filter_annotation.nf` | new: `FILTER_ANNOTATION` |
| `assets/schema_input.json` | `target` ⇒ no `gff3` rule; touch up the two existing `errorMessage`s (§4) |
| `assets/samplesheet_full.csv` | `Lycaena_phlaeas_282391`: `target` + `gff3` → `both`, or blank the `gff3` (§4) |
| `workflows/fomo.nf` | BENCHMARKING gets `PREPROCESSING.out.filtered_gff3` instead of `feature_types`; update the role and feature-type comments |
| `subworkflows/local/preprocessing.nf` | `FILTER_LNC_GFF` → `FILTER_ANNOTATION` (incl. mRNA); `FILTER_TRANSCRIPT.out.mqc` into `ch_mqc_files`; decoy comment; emit comments |
| `subworkflows/local/benchmarking.nf` | remove `FILTER_TARGET` + `feature_types` take; new `ch_annotation` take; refs = role target/both |
| `conf/modules.config` | `ext.args` on `FILTER_TRANSCRIPT`; `FILTER_LNC_GFF` block → `FILTER_ANNOTATION`; delete `FILTER_TARGET`; `sources/` publish rule; publishing header |
| `nextflow.config`, `nextflow_schema.json` | `longest_isoform` (boolean, default `true`) — both, or launch aborts (`failUnrecognisedParams`) |
| `CLAUDE.md` | "Samplesheet & roles" (a `gff3` only on `source`/`both`; pure targets are never benchmarked; replaces "optional for a pure `target`"); results tree (`sources/`); subworkflow map (BENCHMARKING inputs); task counts; Filters; self-pair paragraph; decoy rationale; remove "never TD2-filter the reference"; the id-map contract next to "The merge/split contract" |
| `README.md` | samplesheet rules (`target` ⇒ no `gff3`); outputs: `sources/<species>/annotation/` |
| `bin/run_summary_tables.py` | source funnel: `dropped_non_longest` column + two consistency warnings (§8); the `FILTER_LNC_GFF` comment at :418 |
| `bin/td2_coding_filter.py` | comment only: the `--stage input` branch's `FILTER_LNC_GFF` rationale (its `coding_ids` lose their only consumer; code unchanged) |

No change expected in `projection.nf`, `consensus_top.nf`, `reporting.nf`.
`bin/run_summary_tables.py:418` cross-checks TD2's kept-lncRNA count against the filtered
GFF3 per species — now guaranteed equal by `FILTER_ANNOTATION`; keep it as a regression check.

## Verification

1. **`filter_transcript.sh`, on the five test GFF3s** (`--longest`, `lnc_RNA` and `mRNA`):
   - ≤ 1 transcript per gene; transcript count == gene count;
   - each winner's summed exon length is the gene's max over spliced isoforms — checked by
     an independent pass (`gffread --table @id,@geneid,@covlen` or Python), not the same awk;
   - two runs byte-identical; the tie genes from the profile pick the smallest ID;
   - without `--longest`, byte-identical to `main`'s script.
2. **Fixture GFF3**: a length tie; a gene whose longest isoform is mono-exonic; a
   transcript with no `Parent`; an exon with `Parent=a,b`; a gene with both `mRNA` and
   `lnc_RNA` children; a transcript id with two colons (`transcript:foo:bar`) — the case
   today's normalisation gets wrong.
3. **`filter_gff_by_id.py` keep mode**: each of the three exit-1 checks fires on a crafted
   input; drop mode byte-identical to `main` on a `FILTER_ALLMODELS` input.
4. **Pipeline** (method of plans/18, baseline = `main` via `git stash`):
   - `nextflow lint .` clean on every touched file;
   - `-profile test,singularity --longest_isoform false`: `spliced_fasta`, SeqKit rows and
     the source-side `<sp>.lnc_RNA.filtered.gff3` byte-identical to baseline (proves keep
     mode reproduces today's drop-mode result on real data). gffcompare numbers **will**
     differ (reference now TD2-filtered); record the per-target delta here;
   - `-profile test,singularity` (default): record per-target lncRNA transcript-level
     Sn/Pr/F1 vs. baseline and whether the top-N ranking changed;
   - for each `both` species, the reference staged in `GFFCOMPARE_BATCH`'s work dir has
     the same md5 as `sources/<sp>/annotation/<sp>.lnc_RNA.filtered.gff3`;
   - once with `--include_mrna --include_decoy`: mRNA ★ published, decoys built from
     `FILTER_TRANSCRIPT` output, `decoy_cap` respected;
   - no `targets/null/` or `sources/null/`; `sources/<sp>/` exists for every `source`/`both`
     row and not for the pure target `Lycaena_hippothoe_580924`;
   - **samplesheet rule, both ways**: a local copy of `assets/samplesheet.csv` with a `gff3`
     added to the `target` row aborts at launch with the new message; the committed sheet
     (empty `gff3` on that row) still validates; and `assets/samplesheet_full.csv` validates
     after its fix;
   - source funnel: for every source `dropped_single_exon + dropped_non_longest +
     dropped_coding + kept == raw_lncRNA_transcripts`, `dropped_non_longest` matches the
     multi-isoform profile above (e.g. 5 for `Lycaena_alciphron_282377`: 67 lnc_RNA → 62
     genes), is 0 with `--longest_isoform false`, and neither consistency warning fires;
     the `transcript_filter` section renders once (not twice) in each per-target report.

## Out of scope / possible follow-ups

- Publishing the matching spliced FASTA next to ★ (it is now exactly ★'s sequences).
- Publishing decoy GFF3s (`RELOCATE_LOCI`) under `sources/<sp>/decoy/`.
- Sharing one gunzip between `GUNZIP_FASTA` and `PROJECTION:GUNZIP_TARGET` for `both`
  species (pre-existing duplication).
- A benchmark-only `reference` role, if targets ever need scoring without donating (§4).
- Fixing `RENAME_FASTA_HEADERS`'s `a[2]` colon split (changes every projected id).
- CDS-length ranking for mRNA, or preferring `tag=Ensembl_canonical` before the ID
  tie-break (they already agree everywhere on the lncRNA test data).

## Addendum (2026-10-01): unstranded / trans-spliced transcripts

The first real run after implementation (Apiaceae, `--include_mrna --include_decoy`) died in
`PREPROCESSING:EXTRACT_SEQUENCES` for *Daucus carota* (RefSeq GCF_001625215.2): gffread
`Error parsing strand (?)` on the mitochondrial *nad2* mRNA. **Pre-existing, not caused by
this plan**: `main`'s `filter_transcript.sh` keeps the same three `?`-strand mRNAs (*nad1*,
*nad2*, *rps12*: trans-spliced, exons on both strands); no earlier run had hit a RefSeq
plant annotation with `--include_mrna`.

Fix, in `FILTER_TRANSCRIPT` (the only producer of what gffread extracts and of the
reference): a transcript is kept only if its own strand is `+`/`-`, every child line is on
that same strand, and no parent has an undefined strand. The funnel gains `n_stranded`
(between `n_transcripts` and `n_spliced`) and the run summary a `dropped_unstranded`
column. Verified: test-data outputs byte-identical (with and without `--longest`); on the
real Daucus GFF3 exactly those 3 mRNAs are dropped (49 321 → 49 318), 0 lncRNA affected,
and gffread 0.12.7 (the pipeline's container) then extracts all 27 341 mRNA cleanly.
`GFF_TO_GENE_BED` needs no change: its BED carries no strand, and those loci should stay
masked from decoy placement anyway.

The second error in that run's log — `Join mismatch` from `FILTER_ANNOTATION`'s
`failOnMismatch: true` joins — is a side effect of the abort, not a separate bug: when the
session is torn down the join closes with sources still in flight (their candidates
arrived, their TD2/rename had not), and reports them.

**Resume caveat:** Nextflow does not hash `bin/` script contents into the task key, so
`-resume` reuses the cached (unfixed) `FILTER_TRANSCRIPT` outputs — verified on the test
run. Relaunch an affected run without `-resume` (or delete those task dirs).
