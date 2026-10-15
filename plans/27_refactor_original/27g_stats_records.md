# 27g — Producers emit data records; per-source identity from a manifest

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: U7, T1, T2, T3.
Review items: P1, D3, R1 (structural fix), F2/F7 of the module review.

## Status: not started

## Goal

Every statistics producer also writes **data records** (T1) in the same task. This step is
**additive**: today's `*_mqc.tsv` outputs, adapters and consumers stay untouched. 27h and 27i
switch the consumers and then delete the old path. Per-source identity comes from a manifest
written by the splitter (T3), never from a filename. A bridge test proves the records carry
everything the old TSVs carry.

## Record format (T1)

`<prefix>.<section>.data.json`:

```json
{"schema": "fomo-data/1",
 "section": "gffstats",
 "identity": {"side": "target", "target": "Lycaena_alciphron_282377", "source": "Lycaena_thersamon_265360",
              "track": "lnc_RNA", "label": "lncRNA", "order": 1, "kind": "projected"},
 "data": { …tool-specific payload… }}
```

- A file may hold one record or a JSON list of records.
- Identity keys are drawn from `species, side, target, source, track, label, order, kind,
  stage`. The keys present depend on the section; each section's keys are documented in
  `fomo_stats.py` and tested.
- `fomo_stats.py` gains `make_record`, `write_records`, `read_records(dir)` (validates the
  schema and section) and `SECTIONS` (section → required identity keys).
- Groovy builds the identity JSON once per task in the module script block:
  `groovy.json.JsonOutput.toJson(meta.subMap([...]) + extra)`.
  Batch tasks take per-file identity from the manifest.

## Producers

| section | producer (task) | how the record is written | payload |
|:---|:---|:---|:---|
| `transcript_filter` | FILTER_TRANSCRIPT | `filter_transcript.py --record` (Python, 27f) | n_transcripts, n_stranded, n_spliced, n_selected |
| `gffstats` | GFF_STATS, GFF_STATS_TARGET, GFF_STATS_PROJECTED_BATCH | a bash `printf` wrapper around gff-feature-stats' JSON (works in any image; no parsing) | the gff-feature-stats JSON verbatim |
| `td2_coding` | TD2_CODING_FILTER | Python | n_in, n_coding, n_kept, stage |
| `alignment` | SAMTOOLS_STATS | `bin/samtools_stats_record.awk` (the samtools image has no Python, T2) | the `SN` fields that `samtools_to_mqc.py` reads today, plus reads mapped % |
| `fasta_stats` | SEQKIT_STATS | `bin/seqkit_stats_record.awk` | the seqkit columns `seqkit_to_mqc` keeps today |
| `accuracy` | GFFCOMPARE_BATCH | Python driver (D below) | for each level: Sn, Pr; the query/ref counts; missed/novel counts if present |
| `curation` | CURATE_MODELS | Python (`curate_models.py --record`, or a 5-line wrapper) | `observed`, `discarded_by_reference`, `novel` (curation.json without `params`) |

## Steps

### A. Record helpers and tests

- [ ] **A1. Add the record helpers to `fomo_stats.py`.** pytest covers a round-trip, a missing
  identity key → `FomoError`, and an unknown section → `FomoError`.
- [ ] **A2. Write the two awk emitters** (`bin/samtools_stats_record.awk`,
  `bin/seqkit_stats_record.awk`), with pytest on captured tool outputs (samtools `.stats` and
  seqkit TSV from a V1 run). The emitted JSON must parse, and its numbers must equal the values
  in the old adapters' TSVs for the same input.

### B. Manifest (T3)

- [ ] **B1. `gff_by_source.py --manifest manifest.tsv`** writes `file, source` for every
  per-source file. It writes a header-only manifest when the input is empty.
  `SPLIT_GFF_BY_SOURCE` emits it. pytest covers the manifest rows and the empty case.
- [ ] **B2. Batch inputs carry the manifest.**
  - The aggregate model is an **explicit** manifest row (`file=<target>.allModels.<track>.raw.gff3`,
    `source=allModels_raw`), added where the batch input is built in `projection.nf`.
  - The `ln -s … from_allModels_raw…` symlink trick stays only if a consumer still needs that
    filename. Today the stats filenames are a deliverable (published gffcompare files), so the
    published names must not change (U1).

### C. Emit records (additive)

- [ ] **C1. FILTER_TRANSCRIPT, TD2_CODING_FILTER and CURATE_MODELS** emit `data` (records).
- [ ] **C2. GFF_STATS, GFF_STATS_TARGET and GFF_STATS_PROJECTED_BATCH** emit records with the
  printf wrapper.
  - GFF_STATS_PROJECTED_BATCH takes per-file identity from the manifest and **deletes its bash
    `${id%%.*}` filename parser**.
  - It loops sequentially, under strict mode, so a failing file fails the task with its own exit
    code. Drop `xargs -P`: gff-feature-stats takes about 1 s per file.
- [ ] **C3. SAMTOOLS_STATS and SEQKIT_STATS** emit records with the awk emitters (`img_samtools`,
  `img_seqkit`).

### D. GFFCOMPARE_BATCH

- [ ] **D1. Replace the bash loop with `bin/gffcompare_batch.py`**
  (`--manifest --reference --args --target --track`). It:
  - runs gffcompare sequentially per manifest row with `check=True`;
  - writes **the same published files with the same names** (`<target>.from_<source>.<track>.gffcompare.{stats,tracking,loci,tmap,refmap}`), because these are deliverables;
  - writes one `accuracy` record per comparison using `fomo_stats.parse_accuracy`.
- [ ] **D2. Move the module to `img_curate`** (Python + gffcompare 0.12.6, the same version as
  today's image). G4 must show identical gffcompare deliverables (`#` lines are normalised).
- [ ] **D3. Simplify the per-source re-key in `benchmarking.nf`** (already `transpose()` since
  27b) and keep it working. The stats files still feed REPORTING and RUN_SUMMARY until 27h/27i.

### E. Bridge test and stub matrix

- [ ] **E1. Write `tests/records/test_bridge.py`** (integration, V1 + V2 fixtures).
  - Capture, from a real V2 run, every `*_mqc.tsv` produced *and* every `*.data.json` produced
    for the same items, into `tests/records/fixtures/` (prune to keep it ≤ 5 MB).
  - For every TSV row, the matching record — matched by identity, **not** by parsing the
    Sample label — has the same values: same metrics, numerically equal.
  - This is the oracle that lets 27h and 27i delete the TSV path.
- [ ] **E2. Update the stub matrix (S1–S3)** for the new `data` outputs. Task counts do not change
  in 27g.

## Acceptance

- G1, G2, G3 and G4 identical. Reports and the run summary are unchanged, because nothing
  consumes records yet.
- `grep -rn 'id%%' modules` finds nothing.
- The bridge test passes on V1 and V2.
- Every producer in the table emits records, checked in the S2 stub (outputs exist) and in V2
  (records parse).

## Deviations

| id | step | original plan said | what was done | why (evidence) | impact on later subplans |
|:---|:---|:---|:---|:---|:---|

## Results

*(fill in)*
