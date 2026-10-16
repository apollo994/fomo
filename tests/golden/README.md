# Golden deliverable snapshots

`snapshot.py` reduces a results directory to one line per published deliverable (normalised
sha256); `compare.py` diffs two snapshots. The rules are in plans/27_refactor.md §5.

| file | produced by |
|:---|:---|
| `V1.tsv` … `V4.tsv` | `run_variant.sh V1..V4` on the pre-refactor code (commit `c86521c`), run twice each: identical. The `.loci` rows were recomputed from the same outputs when their rule became `sorted_entries` |
| `V*.info.tsv` | informational: row count and columns of every MultiQC data table |
| `realdata_B0.tsv` | the 23 Lycaenidae `both` species (`realdata_B0.samplesheet.csv`), lncRNA + decoy |

`realdata_B0.tsv` snapshots the results of
`nextflow run <repo> -profile crg --input realdata_B0.samplesheet.csv --include_decoy`, which
reproduces the lncRNA outputs of the 61-species run `thirsty_linnaeus` (session `3098a02f…`,
`-profile crg --include_mrna --include_decoy`) byte for byte. The samplesheet's FASTA/GFF3 paths
are on the CRG cluster.

Gate (G4), from the repository root:

```sh
tests/golden/run_variant.sh V1 /path/to/scratch/V1        # V2, V3, V4 likewise
python3 tests/golden/compare.py tests/golden/V1.tsv /path/to/scratch/V1/snapshot.tsv
```
