# Plan: make the mRNA and decoy tracks optional (`--include_mrna` / `--include_decoy`)

> **Implemented.** Sections below are as-written pre-implementation except where marked
> *(revised during implementation)* — one design point (how the enabled feature-type list
> crosses the subworkflow boundary) had to change after a stub run caught it.

## Context

FOMO's real deliverable is candidate lncRNAs on the target. The mRNA track exists as a
**positive control** (mRNA projects well from close relatives, so its Sn/Pr bounds what the
lncRNA track can achieve) and the decoy track as a **negative control / FDR baseline**
(source loci randomly relocated into the source's own intergenic space, so anything that
projects is a false positive). Both are development instruments, not products — yet both are
unconditionally on, and they dominate runtime. Today, with F=2 feature types and decoys, a
run of N sources costs 2N `FILTER_TRANSCRIPT`, 2N `RELOCATE_LOCI`, N `SAMTOOLS_FAIDX` +
N `GFF_TO_GENE_BED` + N `BEDTOOLS_COMPLEMENT`, **4N minimap2 alignments** + 4N
`BAM_TO_GFF` + 4N `SAMTOOLS_STATS`, 4 `GFFCOMPARE_COMBINE` groups and 4N+4 `GFFCOMPARE`
runs. Only N of those alignments (the real lncRNA ones) contribute to the output GFF3.

`feature_type` is currently minted by two independent hardcoded
`Channel.of('lnc_RNA', 'mRNA')` cross-products — `subworkflows/local/preprocessing.nf:31`
(source side) and `subworkflows/local/benchmarking.nf:20` (target-reference side) — which
must stay in sync because `benchmarking.nf:28-33` pairs projections to references with an
inner `combine(by: 0)` on `meta.feature_type`; a mismatch **silently drops** projections
with no error. `decoy: true` is set in exactly one place, `preprocessing.nf:97`. There are
no `skip_*`/`run_*` params anywhere in the repo and `ext.when` is declared in every local
module but never set.

**Outcome:** default run = lncRNA only (N alignments, 1 combine group). `--include_mrna`
and `--include_decoy` restore today's behaviour exactly. As a side effect the two-place
feature-type duplication is removed: the list is derived once in `workflows/fomo.nf` and
passed to both subworkflows, making the source↔target agreement structural rather than a
comment.

## Decisions

| Question | Decision |
|---|---|
| Param shape | **Two opt-in booleans, `include_mrna = false` and `include_decoy = false`.** Matches the one existing boolean, `include_single_exon = false`; self-documenting in `nextflow run . --help`. Rejected a `--feature_types lnc_RNA,mRNA` CSV list — needs string parsing + a schema `pattern`, and there is no third feature type on the roadmap. |
| Decoy scope when both on | **Decoys for every enabled feature type.** `--include_mrna --include_decoy` reproduces today's 4 tracks bit-for-bit; `--include_decoy` alone gives lncRNA + lncRNA.decoy. Backwards-compatible, and keeps the mRNA null baseline available. |
| Where the enabled list is derived | **Once in `workflows/fomo.nf`, as a plain Groovy `List`, passed into PREPROCESSING and BENCHMARKING as a `take:` input and expanded inside a `flatMap` closure.** One definition, so the inner-join hazard becomes impossible by construction. Rejected a shared helper function file — the repo has no `lib/` or functions-module convention and this needs no new file type. *(revised during implementation — the original plan said `Channel.value([...])` + `.combine()`; see §3.)* |
| How the decoy branch is gated | **A plain `if (params.include_decoy) { … } else { … Channel.empty() }` block in `preprocessing.nf`**, assigning the three decoy-side channels in both arms. Rejected `ext.when` — it would still run/queue the tasks' upstream channel wiring and leaves empty publish dirs; a whole-branch skip is cleaner and matches the pipeline's existing channel-level `filter`/`branch` gating style. |
| Decoy placement vs. `--include_mrna` | **Decoy destinations must stay the complement of the *whole* source annotation — mRNA genes included — even when mRNA is not transferred. Already true; preserved by gating the chain on `include_decoy` alone.** See §4c: the intergenic BED is built from the raw samplesheet GFF3, once per source, and is structurally independent of the feature-type list. |
| Target-reference scope | **BENCHMARKING filters the target annotation only for enabled feature types.** One rule everywhere; with mRNA off the report shows target `raw` + `lncRNA` rows only. |
| `-profile test` | **`conf/test.config` sets both flags true**, so `-profile test` keeps exercising every code path exactly as today and the default path is the lean one. |
| MultiQC row-ordering ladders | **Left untouched.** The `1_raw / 2_lncRNA / 3_decoy_lncRNA / 4_mRNA / 5_decoy_mRNA` ladders in `conf/modules.config` are sparse-safe — a disabled class simply produces no row, and the `'unknown'` fallback can't trigger because the values are always a subset. Renumbering dynamically would churn sample names between runs. |
| `conf/modules.config` | **No changes at all.** Every `meta.decoy` dereference is Groovy-null-safe and the real track keeps `decoy: false`; selectors for skipped processes simply never match. |

## What changes, by cardinality

With N sources, F = enabled feature types (1 or 2), D = 2 if `--include_decoy` else 1:

```
                          today (F=2,D=2)   default (F=1,D=1)   both flags (F=2,D=2)
FILTER_TRANSCRIPT              2N                 N                    2N
RELOCATE_LOCI                  2N                 0                    2N
SAMTOOLS_FAIDX / GENE_BED /     N                 0                     N
  BEDTOOLS_COMPLEMENT
MINIMAP2_ALIGN (+ stats,       4N                 N                    4N
  BAM_TO_GFF)
GFFCOMPARE_COMBINE groups       4                 1                     4
GFFCOMPARE (benchmark)       4N+4               N+1                  4N+4
TD2 lanes                    2 (input+proj)    2 (unchanged)      2 (unchanged)
```

Generalised: `F·N` filters, `F·N·D` alignments, `F·D` combine groups, `F·N·D + F·D`
gffcompares, `F` target references.

## 1. `nextflow.config`

Add the two booleans next to `include_single_exon`:

```groovy
params {
    input            = null
    outdir           = 'results'
    publish_dir_mode = 'copy'
    relocate_seed    = 42
    decoy_cap        = 1000
    include_mrna            = false
    include_decoy           = false
    include_single_exon     = false
    top_consensus_n         = 3
    td2_psauron_min         = 0.5
    ...
}
```

## 2. `nextflow_schema.json`

Mandatory — `validation.failUnrecognisedParams = true` (`nextflow.config:15`) aborts the
run otherwise. Two properties in `$defs.pipeline_options.properties`, placed before
`include_single_exon`, in the existing style (`description` + `help_text` + `fa_icon`):

```json
"include_mrna": {
    "type": "boolean",
    "default": false,
    "description": "Also transfer the source mRNA track (positive control).",
    "help_text": "Off by default: FOMO's deliverable is lncRNA. The mRNA track is a positive control — mRNA projects well from closely related sources, so its gffcompare Sn/Pr bound what the lncRNA track can reach. Turning it on roughly doubles the number of filtering, alignment and benchmarking tasks and adds mRNA rows to every report table.",
    "fa_icon": "fas fa-dna"
},
"include_decoy": {
    "type": "boolean",
    "default": false,
    "description": "Also build and project the relocated decoy track (false-positive baseline).",
    "help_text": "Off by default. Decoys are source loci randomly relocated into that source's own intergenic intervals (see --relocate_seed / --decoy_cap), so any decoy model that projects onto the target is a false positive — the null baseline for the accuracy scatter. One decoy track is built per enabled feature type, so with --include_mrna this doubles again.",
    "fa_icon": "fas fa-mask"
}
```

Also widen the group description at `nextflow_schema.json:53` — "Parameters controlling
decoy generation and the MultiQC report." → "Parameters controlling which feature tracks
are transferred, decoy generation and the MultiQC report." — and append to the existing
`relocate_seed` / `decoy_cap` `help_text` that they take effect only with `--include_decoy`.

## 3. `workflows/fomo.nf` — derive the feature-type list once

Insert after `log.info paramsSummaryLog(workflow)` (`:14`) a plain Groovy `List` (see the
in-file comment for the full rationale):

```groovy
    feature_types = ['lnc_RNA'] + (params.include_mrna ? ['mRNA'] : [])
```

Then thread it through the two call sites:

```groovy
    PREPROCESSING(ch_input.source, feature_types)
    ...
    BENCHMARKING(ch_input.target, PROJECTION.out.gff3, feature_types)
```

**A `List`, not a channel — revised during implementation.** The original plan called for
`Channel.value([...])` combined into each fan-out with `.combine()`. That is **wrong** and
fails silently: `combine` **spreads** a List-valued channel into the tuple, so the
`flatMap { meta, _fa, gff3, ftypes -> ftypes.collect { … } }` closure receives the bare
String `'lnc_RNA'` as `ftypes`, and Groovy's `String.collect{}` iterates its **characters**.
The first stub run produced 7 tracks per source named `l`, `n`, `c`, `_`, `R`, `N`, `A` —
28 alignments instead of 4, no error anywhere. A plain `List` passed to a subworkflow
`take:` stays a `List`, is independently usable by both subworkflows, and needs no
`combine` at all (verified with a standalone two-subworkflow probe, and by the task counts
in Verification 2/3). The `× 7` fan-out is the reason `ch_feature_types` was renamed
`feature_types`: the `ch_` prefix would advertise a channel it is not.

`PROJECTION`, `CONSENSUS_TOP` and `REPORTING` call sites are **unchanged** — they only read
`meta.feature_type`/`meta.decoy` and group by them, so a smaller set of values needs no code
change. The `ch_projected_real` filter at `:36-37` (`!meta.decoy && meta.id != 'combined'`)
stays correct with zero decoys.

## 4. `subworkflows/local/preprocessing.nf`

**4a. New `take:` input and the fan-out (replaces `:29-36`).**

```groovy
    take:
    ch_sources        // [ meta, fasta, gff3 ]
    feature_types     // plain List<String> of enabled feature types (lnc_RNA [, mRNA])
```

```groovy
    // Filter each source GFF3 per enabled feature_type (lnc_RNA [, mRNA])
    ch_filter_input = ch_sources
        .flatMap { meta, _fasta, gff3 ->
            feature_types.collect { ft -> tuple(meta + [feature_type: ft, decoy: false], gff3, ft) }
        }

    FILTER_TRANSCRIPT(ch_filter_input)
```

**4b. Gate the whole decoy branch (replaces `:53-115`).** `SAMTOOLS_FAIDX`,
`GFF_TO_GENE_BED` and `BEDTOOLS_COMPLEMENT` exist *only* to feed `RELOCATE_LOCI` — verified,
`BEDTOOLS_COMPLEMENT.out.bed` is consumed at `:94` and re-emitted as the consumer-less
`intergenic_bed`, and `SAMTOOLS_FAIDX.out.sizes` only at `:66,78` — so all five processes
go inside the guard. Keep the existing channel-construction bodies verbatim; only the
`if`/`else` wrapper and the three channel assignments are new:

```groovy
    // ── Decoy track (opt-in: --include_decoy) ────────────────────────────────
    // Decoys are a false-positive baseline: source loci relocated into that
    // source's own intergenic space. The intergenic machinery (faidx → gene BED
    // → complement) exists only to feed RELOCATE_LOCI, so the whole branch is
    // skipped together. One decoy track per enabled feature_type.
    if (params.include_decoy) {
        SAMTOOLS_FAIDX(...)          // unchanged bodies from :53-115
        GFF_TO_GENE_BED(...)
        BEDTOOLS_COMPLEMENT(...)
        RELOCATE_LOCI(ch_relocate)   // ch_relocate still sets meta + [decoy: true]
        EXTRACT_DECOY_SEQUENCES(ch_decoy_extract.gff, ch_decoy_extract.fasta)

        ch_decoy_gff3      = RELOCATE_LOCI.out.gff3
        ch_decoy_fasta_raw = EXTRACT_DECOY_SEQUENCES.out.gffread_fasta
        ch_intergenic_bed  = BEDTOOLS_COMPLEMENT.out.bed
    }
    else {
        ch_decoy_gff3      = Channel.empty()
        ch_decoy_fasta_raw = Channel.empty()
        ch_intergenic_bed  = Channel.empty()
    }
```

Note `RELOCATE_LOCI` keeps consuming the **pre-TD2** `FILTER_TRANSCRIPT.out.gff3` (`:91`) —
the deliberate unfiltered-null-baseline behaviour documented at `:143-146` is preserved.

**4c. Invariant to preserve — decoy destinations ignore the feature-type list.** Decoys
must land in intergenic space defined against the **entire** source annotation, so a
"decoy" never overlaps a real mRNA gene even when mRNA is not being transferred. This is
already how the code behaves and the guard above keeps it that way; the reason it holds is
worth writing down, because the two obvious "optimisations" both break it:

- `GFF_TO_GENE_BED` is fed from **`ch_sources`** (`preprocessing.nf:63-64`) — the raw
  samplesheet GFF3 — **not** from `FILTER_TRANSCRIPT.out.gff3`. Its awk selector
  (`gff_to_gene_bed.nf:33`) keeps `gene || pseudogene || ncRNA_gene`, which under Ensembl
  GFF3 is protein-coding genes + pseudogenes + all non-coding genes (RefSeq puts lncRNA
  genes under plain `gene`, also covered). `BEDTOOLS_COMPLEMENT` subtracts that merged
  union from the chromosome sizes.
- The chain is **per source, not per feature type** — `tag "${meta.id}"`,
  `ext.prefix = { "${meta.id}" }` (`conf/modules.config:11`), N tasks not F·N — and
  `ch_relocate` (`:91-98`) broadcasts the one BED per species to every feature type via
  `combine(by: 0)` on `meta.id`.

So the feature-type list governs only **which loci are relocated** (`RELOCATE_LOCI`'s
`--feature-type`, from the filtered GFF3), never **where they land**. Two things must
therefore NOT be done while implementing:

1. Do **not** gate any of `SAMTOOLS_FAIDX` / `GFF_TO_GENE_BED` / `BEDTOOLS_COMPLEMENT` on
   `params.include_mrna` — only on `params.include_decoy`. (The guard in §4b does exactly
   this.)
2. Do **not** "simplify" `GFF_TO_GENE_BED`'s input to the filtered per-feature-type GFF3 to
   avoid re-reading the raw file. That would make intergenic space mean "outside lncRNA
   genes", letting lncRNA decoys land inside real mRNA loci — and it would silently change
   meaning the moment `--include_mrna` flips.

**4d. Replace the three downstream references** to the now-conditional processes:

| Line | Was | Becomes |
|---|---|---|
| `:120-123` | `.mix(EXTRACT_DECOY_SEQUENCES.out.gffread_fasta)` | `.mix(ch_decoy_fasta_raw)` |
| `:173` | `.mix(RELOCATE_LOCI.out.gff3.map { … kind: 'decoy' … })` | `.mix(ch_decoy_gff3.map { … kind: 'decoy' … })` |
| `:189` | `decoy_gff3 = RELOCATE_LOCI.out.gff3` | `decoy_gff3 = ch_decoy_gff3` |
| `:188` | `intergenic_bed = BEDTOOLS_COMPLEMENT.out.bed` | `intergenic_bed = ch_intergenic_bed` |

Everything else in the file works unchanged with a subset: the TD2 `branch` at `:129-134`
routes real lncRNA to TD2 and leaves `rest` (possibly empty — `mix` with an empty channel is
a no-op); the `:147-152` `lnc`/`mrna` branch likewise; `ch_renamed_final.filter { meta.decoy }`
at `:190` yields nothing. The `coding_ids` join at `:157` is keyed on `meta.id` only, which
is safe as long as lncRNA remains the **only** class entering TD2 — still true.

**4e. Update the `× 2N` cardinality comments** in the `emit:` block (`:186-191`) to the
`F·N` / `F·N` (decoy: `F·N` or 0) form, and the `× N` on `intergenic_bed` to `N` or 0.

## 5. `subworkflows/local/benchmarking.nf`

New `take:` input and the same `combine`+`flatMap` swap (replaces `:18-25`):

```groovy
    take:
    ch_target          // [ meta(role:'target'), fasta[.gz], gff3[.gz] ] × 1
    ch_projected_gff3  // [ meta(target_id, id, feature_type, decoy), gff3 ]
    ch_feature_types   // value: List<String> — MUST be the same list PREPROCESSING got

    ...
    // Filter target annotation by enabled feature_type — one GFF3 per class.
    ch_target_filter_in = GUNZIP_TARGET_GFF.out.gunzip
        .combine(ch_feature_types)
        .flatMap { meta, gff3, ftypes ->
            ftypes.collect { ft -> tuple(meta + [feature_type: ft, decoy: false], gff3, ft) }
        }
```

Nothing else changes: the projection↔reference `combine(by: 0)` at `:28-33`, the
`multiMap` at `:35-39`, and the target stats union at `:48-50` are all
feature-type-agnostic. Update the `× 16` / `× 20` / `× 2` cardinality comments at
`:10,61-63` to the `F·N·D` form.

## 6. `conf/test.config`

```groovy
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test using single-chromosome subsampled Lycaena assemblies'

    input  = "${projectDir}/assets/samplesheet.csv"
    outdir = 'results'

    // Exercise every track — the test profile is the full-coverage path, while a
    // bare `nextflow run .` is the lean lncRNA-only default.
    include_mrna  = true
    include_decoy = true
}
```

## 7. Docs

**`assets/multiqc/sections.yml`** — description text only, no routing changes:
- `:17` `raw / lncRNA / decoy_lncRNA / mRNA / decoy_mRNA` → append "(the mRNA and decoy
  rows appear only with `--include_mrna` / `--include_decoy`)".
- `:43` and `:64` — same parenthetical after the "(target × source × feature_type ×
  decoy/not)" phrasing.
- `:76` "mRNA and decoy tracks are not filtered." → "…are not filtered (and are absent
  unless `--include_mrna` / `--include_decoy` are set)."

`assets/multiqc/main.yml` needs **no** change — `extra_fn_clean_exts`,
`report_section_order` and `remove_sections` are all track-agnostic.

**`CLAUDE.md`** — four edits:
- New short "**Feature tracks**" note under *Pipeline Architecture*: lncRNA always
  transferred; mRNA (positive control) and decoy (FP baseline) opt-in via `--include_mrna`
  / `--include_decoy`; decoys are built per enabled feature type.
- Subworkflow map table: mark the decoy columns of `PREPROCESSING` (`decoy_spliced_fasta`,
  `decoy_gff3`, `intergenic_bed`) as empty unless `--include_decoy`.
- "**Always-on filters**" paragraph: `FILTER_TRANSCRIPT` stays always-on;
  `RELOCATE_LOCI` moves to "only when `--include_decoy`".
- `meta`-map contract: note that `kind: 'decoy'` and `decoy: true` occur only with
  `--include_decoy`, and that `feature_type` ranges over the enabled subset which is
  derived once in `workflows/fomo.nf`.

## Files to edit

```
nextflow.config                          # 2 new params
nextflow_schema.json                     # 2 new properties + group/help_text wording
workflows/fomo.nf                        # derive ch_feature_types, thread to 2 subworkflows
subworkflows/local/preprocessing.nf      # take:, flatMap fan-out, if(include_decoy) guard
subworkflows/local/benchmarking.nf       # take:, flatMap fan-out
conf/test.config                         # both flags true
assets/multiqc/sections.yml              # description wording
CLAUDE.md                                # 4 doc edits
plans/12_optional_mrna_decoy.md          # this file
```

No changes to `conf/modules.config`, any module under `modules/`, any script under `bin/`,
`assets/schema_input.json`, `subworkflows/local/{projection,consensus_top,reporting}.nf`, or
`assets/multiqc/main.yml`.

## Verification

There is no nf-test/CI in this repo, so verification is a `-stub-run` sweep plus one real run.

### Results actually obtained (4 sources, 1 target, `assets/samplesheet.csv`)

Runs used `-profile singularity` (no docker on this host). **Note the local box has 1 usable
CPU**, so `resourceLimits = [cpus: 1, …]` was needed and any real run is slow.

| Check | Result |
|---|---|
| 1. `--help` renders both params; `--include_mrnaa` aborts | ✅ |
| 2. Lean default stub: 4 alignments (not 16), 1 combine group (not 4), 5 gffcompares, zero `*mRNA*`/`*decoy*` artefacts, and `WARN: no process matching config selector` for `GFF_TO_GENE_BED` / `BEDTOOLS_COMPLEMENT` / `RELOCATE_LOCI` / `EXTRACT_DECOY_SEQUENCES` | ✅ |
| 3. Both flags on vs. `main` (git worktree, same stub config): **310 tasks each, identical name+tag sets, identical 200-file published output set** | ✅ backwards-compatible |
| 4. `--include_decoy` alone: 185 tasks; `--include_mrna false` on the CLI correctly overrides `test.config` (nf-schema coerces the boolean — same 185 tasks as a config-file `false`) | ✅ |
| 6. Decoy placement invariant — see below | ✅ |
| 5, 7. Full real run through MultiQC / lncRNA byte-equivalence | ⏳ **not completed** — the real run is SIGKILLed at ~10 min on this 1-CPU host (twice, including detached via `setsid`). It got as far as `BENCHMARKING:GFFCOMPARE` and was resumable. Run on the cluster (`-profile crg`) to close these two. |

**Verification 6 in detail** (real, non-stub artefacts — the preprocessing stage *did*
complete before the kill, so these are genuine outputs, not stubs):

- `GFF_TO_GENE_BED` ran **4 times = N**, not 8 = F·N, and its staged input GFF3 is a symlink
  straight to `assets/test_data/<sp>/<acc>/<sp>.gff3.gz` — the raw samplesheet file, not any
  `FILTER_TRANSCRIPT` output.
- The pipeline's `Lycaena_hippothoe_580924.genes.bed` from the **mRNA-off** run is
  **byte-identical** (`cmp`, 918 merged intervals) to a BED computed by hand from the full
  raw annotation with the module's own awk. The raw file has 920 `gene` + 133 `ncRNA_gene`
  and no other `*_gene` types, all of which the whitelist covers.
- The resulting intergenic BED (16.2 Mbp) has **zero** overlap with the 920 protein-coding
  genes (13.5 Mbp carved out) and zero with the 133 ncRNA genes.
- Across all 4 sources, **301 real lncRNA decoys**, and with `bedtools intersect`:
  `outside_intergenic = 0` (full containment, `-f 1.0`), `overlap_any_gene = 0`,
  `overlap_CODING_gene = 0` — with mRNA transfer **off**.
- MultiQC row classes emitted by `GFF_STATS_TO_MQC` (what the report renders): `1_raw`,
  `2_lncRNA`, `3_decoy_lncRNA` only — **no `4_mRNA`, no `5_decoy_mRNA`, no `unknown`**,
  confirming the ordering ladders are sparse-safe as designed. (`1_raw.mRNA` rows do appear
  and are correct: those are the *transcript-type* breakdown of the unfiltered input
  annotation, not a transferred track.)

### Commands

1. **Schema wiring.** `nextflow run . --help` — confirm `--include_mrna` / `--include_decoy`
   appear under "Pipeline options". Then `nextflow run . -profile test,docker --include_mrnaa`
   must abort on the unrecognised param (proves `validateParameters` still bites).
2. **Default (lncRNA only), fast structural check.**
   `nextflow run . -profile test,docker -stub-run --include_mrna false --include_decoy false`
   then over the trace `results/pipeline_info/execution_trace.txt`:
   - `grep -c RELOCATE_LOCI` → **0**; likewise `BEDTOOLS_COMPLEMENT`, `GFF_TO_GENE_BED`,
     `SAMTOOLS_FAIDX`, `EXTRACT_DECOY_SEQUENCES` → 0.
   - `grep -c MINIMAP2_ALIGN` → **N** (4 with the test sheet), not 16.
   - `grep -c 'GFFCOMPARE_COMBINE'` → **1**, not 4.
   - no task tag contains `mRNA` or `decoy`.
3. **Both flags on = today's behaviour.** `nextflow run . -profile test,docker -stub-run`
   (test.config turns both on) and diff the sorted process/tag column of its trace against a
   trace captured from `main` before the change — they must be **identical**. This is the
   backwards-compatibility gate.
4. **Decoy-only and mRNA-only combinations** (`-stub-run`, spot-check task counts against the
   cardinality table above):
   - `--include_decoy` → 2N alignments, tags `lnc_RNA` + `lnc_RNA.decoy`, 2 combine groups.
   - `--include_mrna` → 2N alignments, tags `lnc_RNA` + `mRNA`, 2 combine groups, and
     **2** `FILTER_TARGET` tasks; confirm every projection got a gffcompare (N+1 real +
     N+1 mRNA), i.e. nothing was silently dropped by the reference inner join.
5. **Real default run.** `nextflow run . -profile test,docker` and inspect
   `results/multiqc/multiqc_report.html`:
   - `gffstats_input` table: `1_raw` + `2_lncRNA` rows per source, `1_raw` + `2_lncRNA` for
     the target; no `3_decoy_lncRNA` / `4_mRNA` / `5_decoy_mRNA` rows, and **no `unknown`
     row** (that would mean an ordering ladder mis-fired).
   - `gffstats_projection`, `seqkit_stats`, `samtools_align` tables: lncRNA rows only,
     including `from_combined` and `from_top3`.
   - custom accuracy scatter renders with circles only (no ▲ decoy / ■ mRNA markers) and no
     traceback from `bin/gffcompare_accuracy_mqc.py`.
   - `SELECT_TOP_SOURCES` still ranks and emits its table — it filters to real
     non-decoy lncRNA anyway (`bin/select_top_sources.py:78-82`) and errors only when *no*
     real lncRNA stats exist, which can't happen.
6. **Decoy placement is mRNA-aware regardless of `--include_mrna`** (the §4c invariant).
   Run `nextflow run . -profile test,docker --include_mrna false` (decoys on via test.config,
   mRNA off) and check, per source, in the `BEDTOOLS_COMPLEMENT` / `RELOCATE_LOCI` work dirs:
   - The gene BED is built from the **whole** annotation: its interval count must equal the
     merged count of all gene-level features in the raw source GFF3 —
     `zcat <source>.gff3.gz | awk -F'\t' '!/^#/ && ($3=="gene"||$3=="pseudogene"||$3=="ncRNA_gene")' | wc -l`
     should be ≥ `wc -l < <source>.genes.bed` (≥, not =, because overlapping genes are
     merged), and the BED must be **identical** to the one produced by the both-flags-on run
     of step 3 — `cmp` them. Identical BEDs prove the feature-type list cannot reach it.
   - Confirm coding genes really are excluded: the raw source GFF3's protein-coding genes
     must not intersect the intergenic BED —
     `zcat <source>.gff3.gz | awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && $3=="gene" {print $1,$4-1,$5}' | bedtools intersect -a <source>.intergenic.bed -b - -u | wc -l` → **0**.
   - And that the decoys themselves respect it: every relocated model in
     `<source>.lnc_RNA.decoy.gff3` must fall inside the intergenic BED —
     `bedtools intersect -a <(awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && $3=="gene" {print $1,$4-1,$5}' <source>.lnc_RNA.decoy.gff3) -b <source>.intergenic.bed -v | wc -l` → **0**.
   - `grep -c GFF_TO_GENE_BED` in the trace → **N**, not F·N (one BED per source, shared
     across feature types).
7. **Output equivalence.** The default run's `results/gffcompare/*.lnc_RNA.*` and
   `results/combined_gtf_to_gff/*` lncRNA products should be byte-identical to the
   corresponding files from a both-flags-on run — the lncRNA track must not be perturbed by
   the presence or absence of its neighbours. (Decoy relocation consumes
   `params.relocate_seed` per-source independently, so no RNG coupling exists; this check
   confirms it.)

## Risks / follow-ups

- **Value-channel `combine` broadcast.** `Channel.value([...])` combined with a queue channel
  must broadcast, not consume-once. This is the same idiom already in use at
  `projection.nf:32-34,40` (`.first()` → value channel → `.combine`), so it is proven in this
  codebase — but verify in step 2 that `FILTER_TRANSCRIPT` runs N times, not once.
- **Nothing enforces the two subworkflows share the list** at the type level; the guard is
  that `workflows/fomo.nf` is the only place that constructs it. Worth a comment on both
  `take:` declarations (included above) so a future edit doesn't re-hardcode one side.
- **Pre-existing, not introduced here: `GFF_TO_GENE_BED`'s gene-type whitelist is
  closed-set.** `gff_to_gene_bed.nf:33` matches exactly `gene`, `pseudogene`, `ncRNA_gene`.
  That covers Ensembl and RefSeq (the two sources in `assets/samplesheet*.csv`), but an
  annotation using other SO gene terms would leave those loci looking "intergenic" and
  decoys could land on them. Step 6's coding-gene intersect check is the detector. If a
  future source needs it, widen the selector to `$3 ~ /_gene$/ || $3=="gene" || $3=="pseudogene"`
  — a separate change, and one that only matters when `--include_decoy` is set.
- **`ext.when` stays unused.** Every local module still carries
  `when: task.ext.when == null || task.ext.when`; this plan does not start using it. If a
  finer-grained per-process skip is ever wanted, that hook is there.
- **Multi-target branch (`multitarget`) will conflict** in `workflows/fomo.nf` and both
  subworkflow `take:` blocks. Small, but rebase deliberately.
- **Not addressed:** the five hardcoded `lnc_RNA` `ext.prefix` closures on the TD2 lanes
  (`conf/modules.config:41,45,168,177,181`) and the `(lnc_RNA|mRNA)` filename regexes in
  `bin/select_top_sources.py:24` / `bin/gffcompare_accuracy_mqc.py:46`. They are correct for
  any subset of the current two classes, so they don't block this change — but they are what
  would need touching to add a *third* feature type.
