# FOMO Pipeline — End-to-End Implementation Plan

## Context

The FOMO pipeline currently takes pre-extracted mRNA FASTA files and runs BUSCO on them (2 steps). The goal is to expand it into an end-to-end pipeline that starts from **genome assemblies + GFF3 annotations** (as in the existing `assets/samplesheet.tsv`), extracts mRNA features, gets sequences, and runs BUSCO — with nf-test coverage for every step.

The architecture follows patterns from `nf-core-crossmap`: proper directory layout, `conf/` configs, nf-core modules, and nf-test for each process.

---

## Pipeline Steps (End-to-End)

```
samplesheet.csv (species, assembly, annotation)
  │
  ├─[1] Parse samplesheet → [meta, assembly, annotation] channels
  │
  ├─[2] FILTER_FEATURES (local module)
  │     → filter GFF3 to keep only mRNA features
  │
  ├─[3] GFFREAD (nf-core module)
  │     → extract mRNA transcript sequences (filtered GFF3 + genome FASTA → mRNA FASTA)
  │     → uses gffread for all sequence extraction (fast, lightweight)
  │
  ├─[4] BUSCO_DOWNLOAD (existing, conditional)
  │     → download lineage DB if not pre-cached
  │
  └─[5] BUSCO_BUSCO (existing)
        → assess transcriptome completeness
        → emit short_summaries_json
```

> **Note:** GFFREAD handles sequence extraction directly from filtered GFF3 — no separate
> longest-isoform step needed. AGAT is only used for feature filtering (step 2).

## Implementation Steps

### Step 1: Project Restructure (following nf-core-crossmap patterns)

**New directory layout:**
```
main.nf                          # Entry point — calls FOMO workflow
workflows/
  fomo.nf                        # Main workflow logic (moved from main.nf)
conf/
  base.config                    # Resource labels (process_single, _low, _medium)
  modules.config                 # Per-module ext.args and publishDir
  test.config                    # Test profile with small data
modules/
  nf-core/
    busco/busco/                 # (existing)
    agat/spkeeplongestisoform/   # (install via nf-core)
    gffread/                     # (install via nf-core)
  local/
    busco_download/              # (existing)
    filter_features/             # Local: extract mRNA from GFF3
tests/
  workflows/
    fomo.nf.test                 # Workflow-level nf-test
```

### Step 2: Install nf-core Modules

```bash
conda activate nf-core
cd /users/rg/fzanarello/pipelines/fomo
nf-core modules install agat/spkeeplongestisoform
nf-core modules install gffread
```

These provide:
- **agat/spkeeplongestisoform** — keeps longest isoform per gene from GFF3
- **gffread** — extracts transcript sequences from GFF3 + genome FASTA

### Step 3: Create Local FILTER_FEATURES Module

**File:** `modules/local/filter_features/main.nf`

Extracts mRNA features from a GFF3 file. Inspired by `nf-core-crossmap/modules/local/filter_features.nf`:
- Uses AGAT container (same as other AGAT modules)
- Input: `[meta, gff3]` + `val feature_type` (default: `mRNA`)
- Output: `[meta, filtered_gff3]`
- Logic: awk to extract mRNA IDs, then `agat_sp_filter_feature_from_keep_list.pl`

**Test:** `modules/local/filter_features/tests/main.nf.test`
- Uses one of the Papilio test GFF3 files
- Asserts process success
- Validates output GFF3 contains only mRNA features

### Step 4: Create conf/ Configuration Files

**`conf/base.config`** — resource labels:
```groovy
process {
    withLabel:process_single { cpus=1; memory=6.GB; time=4.h }
    withLabel:process_low    { cpus=2; memory=12.GB; time=4.h }
    withLabel:process_medium { cpus=6; memory=36.GB; time=8.h }
    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries = 1
}
```

**`conf/modules.config`** — per-module ext.args and publishDir:
```groovy
process {
    withName: 'GFFREAD' { ext.args = '-w'; publishDir = [...] }
    withName: 'BUSCO_BUSCO' { publishDir = [...] }
    // etc.
}
```

**`conf/test.config`** — minimal test profile pointing to `assets/samplesheet.tsv` with Papilio data, using `lepidoptera_odb12` lineage.

### Step 5: Refactor main.nf and Create workflows/fomo.nf

**`main.nf`** (simplified entry point):
```groovy
include { FOMO } from './workflows/fomo'
workflow { FOMO() }
```

**`workflows/fomo.nf`** (full pipeline logic):
1. Parse `samplesheet.tsv` → create `[meta, assembly, annotation]` channel
2. FILTER_FEATURES(annotation, 'mRNA') → filtered GFF3
3. AGAT_SPKEEPLONGESTISOFORM(filtered_gff3) → longest isoform GFF3
4. GFFREAD(longest_gff3, assembly) → mRNA FASTA
5. BUSCO_DOWNLOAD(lineage) — conditional
6. BUSCO_BUSCO(mRNA_fasta, mode, lineage, lineage_path, [], true)

### Step 6: Update nextflow.config

- Add `includeConfig 'conf/base.config'` and `conf/modules.config`
- Change `params.input` to point to samplesheet TSV instead of glob
- Add `test` profile loading `conf/test.config`
- Keep existing profiles (docker, conda, crg)

### Step 7: Write nf-test for Each Step

All tests use the Papilio test data in `assets/test_data/`.

| Test File | What It Tests | Key Assertions |
|-----------|---------------|----------------|
| `modules/local/filter_features/tests/main.nf.test` | mRNA extraction from GFF3 | Success; output has mRNA features only |
| `modules/nf-core/agat/spkeeplongestisoform/tests/` | (ships with nf-core module) | Already included |
| `modules/nf-core/gffread/tests/` | (ships with nf-core module) | Already included |
| `modules/local/busco_download/tests/main.nf.test` | Lineage download | Success; `busco_downloads/` exists |
| `modules/nf-core/busco/busco/tests/` | (existing — 7 tests) | Already included |
| `tests/workflows/fomo.nf.test` | Full end-to-end workflow | All steps succeed; BUSCO JSON output produced |

**Custom tests to write (3 files):**

1. **FILTER_FEATURES test** — feed a test GFF3, assert only mRNA in output
2. **BUSCO_DOWNLOAD test** — run with `bacteria_odb10` (small), assert download dir exists
3. **Workflow test** — run the full FOMO workflow with test profile, assert final BUSCO output

### Step 8: Create nf-test Config

**`nf-test.config`** (root):
```groovy
config {
    testsDir "."
    workDir ".nf-test"
    configFile "nextflow.config"
    profile "docker,test"
}
```

---

## Files to Create/Modify

| Action | File | Description |
|--------|------|-------------|
| **Create** | `workflows/fomo.nf` | Main workflow logic |
| **Create** | `modules/local/filter_features/main.nf` | mRNA filter module |
| **Create** | `modules/local/filter_features/tests/main.nf.test` | nf-test for filter |
| **Create** | `modules/local/filter_features/tests/nextflow.config` | Test config |
| **Create** | `modules/local/busco_download/tests/main.nf.test` | nf-test for download |
| **Create** | `modules/local/busco_download/tests/nextflow.config` | Test config |
| **Create** | `tests/workflows/fomo.nf.test` | Workflow-level test |
| **Create** | `tests/workflows/nextflow.config` | Workflow test config |
| **Create** | `conf/base.config` | Resource labels |
| **Create** | `conf/modules.config` | Module ext.args/publishDir |
| **Create** | `conf/test.config` | Test profile |
| **Create** | `nf-test.config` | nf-test framework config |
| **Modify** | `main.nf` | Simplify to just call workflow |
| **Modify** | `nextflow.config` | Add conf includes, test profile, samplesheet input |
| **Install** | `modules/nf-core/agat/spkeeplongestisoform/` | via `nf-core modules install` |
| **Install** | `modules/nf-core/gffread/` | via `nf-core modules install` |
| **Update** | `CLAUDE.md` | Reflect new architecture |

---

## Verification

```bash
# Activate nf-core environment
conda activate nf-core

# 1. Run individual module tests
nf-test test modules/local/filter_features/tests/main.nf.test --profile docker
nf-test test modules/local/busco_download/tests/main.nf.test --profile docker

# 2. Run nf-core module tests (already ship with the modules)
nf-test test modules/nf-core/agat/spkeeplongestisoform/tests/main.nf.test --profile docker
nf-test test modules/nf-core/gffread/tests/main.nf.test --profile docker
nf-test test modules/nf-core/busco/busco/tests/main.nf.test --profile docker

# 3. Run workflow-level test
nf-test test tests/workflows/fomo.nf.test --profile docker,test

# 4. Run full pipeline manually with test data
nextflow run main.nf -profile docker,test --outdir results
```

---

## Decisions (confirmed with user)

1. **Feature type**: mRNA only — keep the pipeline focused on annotation completeness assessment
2. **Samplesheet format**: Switch to CSV (nf-core convention) — update `get_test_data.sh` accordingly
3. **Filter approach**: Local FILTER_FEATURES module using awk + AGAT keep-list (like crossmap) for full control
