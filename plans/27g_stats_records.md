# 27g — Producers emit data records; per-source identity from a manifest

Part of [27_refactor.md](27_refactor.md) (read §3–§8 there first). Decisions: U7, T1, T2, T3.
Review items: P1, D3, R1 (structural fix), F2/F7 of the module review.

## Status: done (2026-10-15)

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
- (updated by D27g-5: one `ext.identity` closure at the end of `conf/modules.config` builds it for
  every process instead.) Groovy builds the identity JSON once per task in the module script block:
  `groovy.json.JsonOutput.toJson(meta.subMap([...]) + extra)`.
  Batch tasks take per-file identity from the manifest.

## Producers

| section | producer (task) | how the record is written | payload |
|:---|:---|:---|:---|
| `transcript_filter` | FILTER_TRANSCRIPT | `filter_transcript.py --record` (Python, 27f) | n_transcripts, n_stranded, n_spliced, n_selected |
| `gffstats` | GFF_STATS, GFF_STATS_TARGET, GFF_STATS_PROJECTED_BATCH | a bash `printf` wrapper around gff-feature-stats' JSON (works in any image; no parsing) | the gff-feature-stats JSON verbatim |
| `td2_coding` | TD2_CODING_FILTER | Python | n_in, n_coding, n_kept, stage |
| `alignment` | SAMTOOLS_STATS | `bin/samtools_stats_record.awk` (the samtools image has no Python, T2; `.sh`, D27g-1) | the `SN` fields that `samtools_to_mqc.py` reads today, plus reads mapped % |
| `fasta_stats` | SEQKIT_STATS | `bin/seqkit_stats_record.awk` (`.sh`, D27g-1) | the seqkit columns `seqkit_to_mqc` keeps today |
| `accuracy` | GFFCOMPARE_BATCH | Python driver (D below; a bash driver + `bin/gffcompare_stats_record.sh`, D27g-3) | for each level: Sn, Pr; the query/ref counts; missed/novel counts if present |
| `curation` | CURATE_MODELS | Python (`curate_models.py --record`, or a 5-line wrapper) | `observed`, `discarded_by_reference`, `novel` (curation.json without `params`) |

## Steps

### A. Record helpers and tests

- [x] **A1. Add the record helpers to `fomo_stats.py`.** pytest covers a round-trip, a missing
  identity key → `FomoError`, and an unknown section → `FomoError`.
- [x] **A2. Write the two awk emitters** (`bin/samtools_stats_record.sh`,
  `bin/seqkit_stats_record.sh`; updated by D27g-1), with pytest on captured tool outputs (samtools `.stats` and
  seqkit TSV from a V1 run). The emitted JSON must parse, and its numbers must equal the values
  in the old adapters' TSVs for the same input.

### B. Manifest (T3)

- [x] **B1. `gff_by_source.py --manifest manifest.tsv`** writes `file, source` for every
  per-source file. It writes a header-only manifest when the input is empty.
  `SPLIT_GFF_BY_SOURCE` emits it. pytest covers the manifest rows and the empty case.
- [x] **B2. Batch inputs carry the manifest.** (updated by D27g-2)
  - The aggregate model is an **explicit** manifest row (`file=<target>.allModels.<track>.raw.gff3`,
    `source=allModels_raw`), added where the batch input is built in `projection.nf`.
  - The `ln -s … from_allModels_raw…` symlink trick stays only if a consumer still needs that
    filename. Today the stats filenames are a deliverable (published gffcompare files), so the
    published names must not change (U1).

### C. Emit records (additive)

- [x] **C1. FILTER_TRANSCRIPT, TD2_CODING_FILTER and CURATE_MODELS** emit `data` (records).
- [x] **C2. GFF_STATS, GFF_STATS_TARGET and GFF_STATS_PROJECTED_BATCH** (updated by D27g-4) emit records with the
  printf wrapper.
  - GFF_STATS_PROJECTED_BATCH takes per-file identity from the manifest and **deletes its bash
    `${id%%.*}` filename parser**.
  - It loops sequentially, under strict mode, so a failing file fails the task with its own exit
    code. Drop `xargs -P`: gff-feature-stats takes about 1 s per file.
- [x] **C3. SAMTOOLS_STATS and SEQKIT_STATS** emit records with the awk emitters (`img_samtools`,
  `img_seqkit`).

### D. GFFCOMPARE_BATCH

- [x] **D1. Replace the bash loop with `bin/gffcompare_batch.py`** (updated by D27g-3: `bin/gffcompare_batch.sh` + `bin/gffcompare_stats_record.sh`)
  (`--manifest --reference --args --target --track`). It:
  - runs gffcompare sequentially per manifest row with `check=True`;
  - writes **the same published files with the same names** (`<target>.from_<source>.<track>.gffcompare.{stats,tracking,loci,tmap,refmap}`), because these are deliverables;
  - writes one `accuracy` record per comparison using `fomo_stats.parse_accuracy`.
- [x] ~~**D2. Move the module to `img_curate`**~~ (not done, D27g-3) (Python + gffcompare 0.12.6, the same version as
  today's image). G4 must show identical gffcompare deliverables (`#` lines are normalised).
- [x] **D3. Simplify the per-source re-key in `benchmarking.nf`** (already `transpose()` since
  27b) and keep it working. The stats files still feed REPORTING and RUN_SUMMARY until 27h/27i.

### E. Bridge test and stub matrix

- [x] **E1. Write `tests/records/test_bridge.py`** (integration, V1 + V2 fixtures).
  - Capture, from a real V2 run, every `*_mqc.tsv` produced *and* every `*.data.json` produced
    for the same items, into `tests/records/fixtures/` (prune to keep it ≤ 5 MB).
  - For every TSV row, the matching record — matched by identity, **not** by parsing the
    Sample label — has the same values: same metrics, numerically equal.
  - This is the oracle that lets 27h and 27i delete the TSV path.
- [x] **E2. Update the stub matrix (S1–S3)** for the new `data` outputs. Task counts do not change
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
| D27g-1 | A2, D1 | awk emitters `bin/*.awk` | bash scripts running awk: `samtools_stats_record.sh`, `seqkit_stats_record.sh`, `gffcompare_stats_record.sh` | as D27c-1 (no portable awk shebang). They avoid the `var (expr)` concatenation, which busybox awk 1.32 (gffcompare image) parses as a function call ("Call to undefined function"); checked under gawk, mawk and both busybox versions | none |
| D27g-2 | B2 | the aggregate manifest row is added where the batch input is built in `projection.nf` | `gff_by_source.py --aggregate-id allModels_raw` writes it as the manifest's first row | the splitter knows the input file; adding a row channel-side means writing a file in a closure. Still one explicit row, not a filename trick | none |
| D27g-3 | D1, D2 | Python driver on `img_curate` | GFFCOMPARE_BATCH stays on `img_gffcompare` (BioContainers 0.12.6); the driver is `bin/gffcompare_batch.sh`, the records come from `bin/gffcompare_stats_record.sh` (equal to `fomo_stats.parse_gffcompare_stats`, tested) | with `img_curate`'s gffcompare (bioconda build h9948957_4) G4 V2 changed 3 published `.loci` files: the same transcripts, tied entries in another order. Run standalone, both builds change that order with the `-o` prefix and the directory (BioContainers: short prefix `e15fc…` = the pipeline's, long prefix `c4517…`), so the order is reproducible only for one build and path layout; G4's gate directories have equal path lengths | 27k G5 runs in another directory than B0: expect `.loci` tie-order differences there and compare them order-insensitively only with evidence (updated by D27g-3) |
| D27g-4 | C2 | drop `xargs -P` | both batch modules loop sequentially; `process_batch` (4 CPUs, retry on 123/125) is deleted and they use `process_low` | a batch task takes 1–3 s in the 23-species run; with one worker the xargs exit codes cannot occur | effective resources of GFFCOMPARE_BATCH / GFF_STATS_PROJECTED_BATCH change (4 → 1 CPU, 2 h → 15 min × attempt) |
| D27g-5 | T1 | Groovy builds the identity in each module | one `ext.identity` closure for every process, at the end of `conf/modules.config` | one definition instead of eight. Declared before a `withName` block, a process-level `ext.*` captures that block's `ext.*` (`nextflow config` warned "Unknown config attribute process.withName:…", nf-test failed); a test pins it after the last `withName`. A user `-c` file is not affected (checked) | 27j documents the ordering rule |

## Results

**Record format.** `fomo_stats.py`: `SCHEMA`, `IDENTITY_KEYS`, `SECTIONS` (section → required
identity keys), `make_record`, `write_records`, `read_records`, `identity_arg`, and
`parse_gffcompare_stats` (levels + query/reference/matching/missed/novel counts). Identity:
`species` (= `meta.id`), `side` (`target` when the task has a `target_id`, else `source`),
`target`, `track`, `label`, `order`, `kind`, plus `source` (batch rows, from the manifest) and
`stage` (TD2, added by the script). The curation record also carries `has_reference`, which
the discarded table depends on.

**Producers** (all additive; the `*_mqc.tsv` path is untouched):

| section | producer | record writer |
|:---|:---|:---|
| `transcript_filter` | FILTER_TRANSCRIPT | `filter_transcript.py --record --identity` |
| `gffstats` | GFF_STATS, GFF_STATS_TARGET, GFF_STATS_PROJECTED_BATCH | printf around the gff-feature-stats JSON |
| `td2_coding` | TD2_CODING_FILTER | `td2_coding_filter.py --identity` |
| `alignment` | SAMTOOLS_STATS | `samtools_stats_record.sh` |
| `fasta_stats` | SEQKIT_STATS | `seqkit_stats_record.sh` |
| `accuracy` | GFFCOMPARE_BATCH | `gffcompare_stats_record.sh` (via `gffcompare_batch.sh`) |
| `curation` | CURATE_MODELS | `curate_models.py --identity` |

Records, manifests and MultiQC tables are kept out of the published directories (saveAs), so
the published path set is unchanged. The four subworkflows emit `data` for 27h.

**Manifest.** SPLIT_GFF_BY_SOURCE emits `<target>.<track>.manifest.tsv` (`file, source`; the
aggregate first). Both batch tasks read each file's source from it; the bash `${id%%.*}`
filename parser is gone (`grep -rn 'id%%' modules` is empty). GFFCOMPARE_BATCH still links the
aggregate as `<target>.from_allModels_raw.<track>.projected.gff3`, because its published
`.tmap`/`.refmap` names embed the query name.

**Bridge (E1).** `tests/records/capture_bridge.py` pairs every MultiQC table with the record of
the same measurement through task lineage; fixtures from the 27g gate runs: V1 33 pairs, V2 88
pairs, all 9 table kinds, 56 KB. `tests/records/test_bridge.py`: every pair holds the same
values (gffstats rows rebuilt with the adapter's own functions, other tables numerically,
accuracy levels and the filename identity equal to the record's).

**V2 records:** 74 files, 74 valid records (gffstats 34, accuracy 8, curation 8, fasta_stats 8,
alignment 8, td2_coding 4, transcript_filter 4).

**Gates:** G1 278 passed (121 of them bridge pairs); G2 pass; G3 S1–S4 pass (task counts unchanged); G4 V1–V4 identical
(after D27g-3). Cluster use 1.9 CPU-h (three gate runs, three stub runs).
