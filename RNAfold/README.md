# RNAfold sliding-window MFE test

Standalone, exploratory experiment — **not** part of the FOMO Nextflow pipeline.

Question: does ViennaRNA secondary-structure stability (minimum free energy, MFE)
distinguish **lncRNA** from **mRNA**, and **real** candidates from relocated **decoys**?

It slides windows of **30 / 60 / 90 nt** (step **10 nt**) over each transcript, computes
the per-window MFE with the ViennaRNA Python API (`RNA.fold_compound(...).mfe()`, see
`helpers.py`), and compares the resulting distributions across classes.

## Setup

```bash
conda env create -f RNAfold/environment.yml
conda activate rnafold_test
```

## Run

From the repository root:

```bash
python RNAfold/run_rnafold_test.py
```

## What it does

- Species: `Lycaena_alciphron_282377` (single species).
- Samples **50 sequences** (fixed seed `42`) from each of 4 classes, drawn from
  `results/rename_fasta_headers/`:
  - lncRNA real / lncRNA decoy / mRNA real / mRNA decoy.
- Writes the sampled subsets to `subsets/<class>_<source>.50.fasta`.

## Outputs

| Path                                   | Contents                                                                 |
|----------------------------------------|--------------------------------------------------------------------------|
| `results/mfe_long.csv`                 | Tidy table: `seq_id, rna_class, source, window_size, center, mfe, mfe_per_nt` |
| `results/mfe_summary.csv`              | Mean/median/std MFE & MFE-per-nt grouped by class × source × window size |
| `plots/mfe_distribution.png`           | Box plots of MFE by class × source, faceted by window size               |
| `plots/mfe_per_nt_distribution.png`    | Same, length-normalised (MFE / window size) for cross-window comparison  |
| `plots/positional_profile.png`         | Mean MFE vs relative position along transcript (5'→3'), per class × source |

## Parameters

Edit the constants at the top of `run_rnafold_test.py` to change them:
`SEED`, `N_PER_CLASS`, `WINDOW_SIZES`, `STEP`, `SPECIES`.
