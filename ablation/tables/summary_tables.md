# Ablation summary tables

## A1 support threshold — lnc_RNA

| k (min species) | supported chains (sum) | Sn strict (median) | Pr strict (median) | F1 strict (median) | Pr lenient (pooled) | decoy chains (sum) | curated after ref filter (sum) |
|:---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 288,410 | 15.7 | 3.7 | 5.2 | 3.8 | 120,504 | 211,121 |
| 2 | 38,934 | 8.4 | 10.1 | 8.5 | 14.7 | 154 | 21,589 |
| 3 | 14,257 | 5.5 | 16.8 | 7.8 | 20.7 | 23 | 6,976 |
| 4 | 5,712 | 4.3 | 29.5 | 7.5 | 33.0 | 5 | 2,093 |
| 5 | 3,006 | 2.1 | 35.7 | 3.9 | 40.7 | 0 | 853 |
| 6 | 1,922 | 1.5 | 42.6 | 2.9 | 47.5 | 0 | 457 |
| 7 | 1,325 | 1.1 | 46.7 | 2.1 | 50.2 | 0 | 273 |
| 8 | 915 | 1.0 | 51.9 | 1.9 | 55.6 | 0 | 153 |
| 9 | 599 | 0.8 | 60.6 | 1.5 | 63.1 | 0 | 66 |
| 10 | 415 | 0.6 | 65.2 | 1.2 | 66.7 | 0 | 36 |
| 11 | 287 | 0.4 | 63.6 | 0.8 | 65.9 | 0 | 22 |
| 12 | 200 | 0.2 | 60.0 | 0.4 | 65.0 | 0 | 14 |

## A1 support threshold — mRNA

| k (min species) | supported chains (sum) | Sn strict (median) | Pr strict (median) | F1 strict (median) | Pr lenient (pooled) | decoy chains (sum) | curated after ref filter (sum) |
|:---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2,400,643 | 60.0 | 7.0 | 12.6 | 8.3 | 242,797 | 164,181 |
| 2 | 612,911 | 51.3 | 24.2 | 32.8 | 25.9 | 405 | 24,212 |
| 3 | 355,120 | 43.3 | 35.6 | 39.3 | 37.9 | 87 | 10,231 |
| 4 | 240,903 | 37.9 | 42.7 | 39.1 | 47.8 | 33 | 5,199 |
| 5 | 179,456 | 33.2 | 50.3 | 39.5 | 54.2 | 17 | 3,112 |
| 6 | 141,760 | 29.8 | 56.2 | 38.5 | 59.8 | 12 | 1,978 |
| 7 | 114,406 | 26.9 | 61.3 | 36.8 | 64.1 | 7 | 1,273 |
| 8 | 93,805 | 23.7 | 66.1 | 34.5 | 68.0 | 4 | 827 |
| 9 | 76,981 | 20.4 | 69.6 | 31.5 | 70.8 | 1 | 564 |
| 10 | 62,595 | 16.9 | 72.9 | 27.7 | 73.4 | 1 | 351 |
| 11 | 50,223 | 13.7 | 74.6 | 23.3 | 75.7 | 0 | 199 |
| 12 | 39,614 | 10.9 | 77.0 | 19.1 | 77.9 | 0 | 99 |

## A2 reference filter (lncRNA, k = 2, sums over targets)

| supported chains | kept (novel) | lost to reference | lost: = annotated lncRNA isoform | lost: = annotated mRNA isoform | lost: touches mRNA exon |
|:---|---:|---:|---:|---:|---:|
| 38,934 | 21,589 | 17,345 | 5,710 | 977 | 7,636 |

## A3 calibration by exact support — lnc_RNA

| species supporting the chain | chains | % = strict reference | % = any annotated isoform | decoy chains |
|:---|---:|---:|---:|---:|
| 1 | 249,476 | 1.9 | 2.1 | 120,350 |
| 2 | 24,677 | 10.6 | 11.2 | 131 |
| 3 | 8,545 | 12.1 | 12.6 | 18 |
| 4 | 2,706 | 23.2 | 24.4 | 5 |
| 5–6 | 1,681 | 30.6 | 33.2 | 0 |
| 7–9 | 910 | 41.0 | 42.6 | 0 |
| 10–14 | 376 | 63.6 | 67.8 | 0 |
| 15–+ | 39 | 46.2 | 56.4 | 0 |

## A3 calibration by exact support — mRNA

| species supporting the chain | chains | % = strict reference | % = any annotated isoform | decoy chains |
|:---|---:|---:|---:|---:|
| 1 | 1,787,732 | 1.6 | 2.2 | 242,392 |
| 2 | 257,791 | 7.4 | 9.5 | 318 |
| 3 | 114,217 | 14.1 | 17.0 | 54 |
| 4 | 61,447 | 25.6 | 29.0 | 16 |
| 5–6 | 65,050 | 32.8 | 36.8 | 10 |
| 7–9 | 51,811 | 48.5 | 52.9 | 6 |
| 10–14 | 44,298 | 65.8 | 69.5 | 1 |
| 15–+ | 18,297 | 80.5 | 82.8 | 0 |

## A4 donor sets — targets with >= 2 congeneric donors (n = 16)

| donor set | targets | chains (median) | Sn strict | Pr strict | F1 strict | F1 lenient |
|:---|---:|---:|---:|---:|---:|---:|
| all 22 donors, k = 1 (union) | 16 | 17,596.0 | 36.5 | 4.5 | 7.7 | 8.0 |
| all 22 donors, k = 2 (pipeline default) | 16 | 2,217.0 | 18.1 | 17.9 | 17.9 | 17.7 |
| nearest donor only, k = 1 | 16 | 1,321.5 | 23.0 | 30.3 | 26.1 | 25.4 |
| nearest 2, k = 2 | 16 | 408.5 | 14.4 | 59.9 | 23.1 | 21.6 |
| nearest 3, k = 2 | 16 | 496.0 | 14.9 | 54.4 | 23.7 | 22.1 |
| nearest 4, k = 2 | 16 | 626.0 | 15.5 | 50.5 | 24.0 | 22.5 |
| nearest 6, k = 2 | 16 | 928.0 | 16.5 | 35.0 | 22.3 | 21.3 |
| nearest 8, k = 2 | 16 | 1,362.0 | 17.1 | 23.8 | 19.8 | 19.2 |
| nearest 12, k = 2 | 16 | 2,002.0 | 17.9 | 20.5 | 18.5 | 18.2 |
| congeners only, k = 2 | 16 | 555.0 | 14.4 | 50.5 | 23.1 | 21.6 |
| no congeners, k = 2 | 16 | 1,513.5 | 5.5 | 7.1 | 5.7 | 5.7 |
| farthest 12, k = 2 | 16 | 325.5 | 1.6 | 13.0 | 2.9 | 3.0 |

## A4 donor sets — targets with < 2 congeneric donors (n = 7)

| donor set | targets | chains (median) | Sn strict | Pr strict | F1 strict | F1 lenient |
|:---|---:|---:|---:|---:|---:|---:|
| all 22 donors, k = 1 (union) | 7 | 9,868.0 | 6.6 | 2.1 | 3.2 | 3.2 |
| all 22 donors, k = 2 (pipeline default) | 7 | 1,315.0 | 3.2 | 5.8 | 4.1 | 4.2 |
| nearest donor only, k = 1 | 7 | 1,503.0 | 4.7 | 6.7 | 5.7 | 5.6 |
| nearest 2, k = 2 | 7 | 40.0 | 0.7 | 44.1 | 1.3 | 1.2 |
| nearest 3, k = 2 | 7 | 75.0 | 1.2 | 29.9 | 2.4 | 2.2 |
| nearest 4, k = 2 | 7 | 139.0 | 1.5 | 20.9 | 2.8 | 2.6 |
| nearest 6, k = 2 | 7 | 403.0 | 2.0 | 12.7 | 3.4 | 3.5 |
| nearest 8, k = 2 | 7 | 471.0 | 2.4 | 9.9 | 3.9 | 4.0 |
| nearest 12, k = 2 | 7 | 805.0 | 2.9 | 7.6 | 4.2 | 4.2 |
| no congeners, k = 2 | 6 | 1,280.5 | 2.5 | 5.0 | 3.2 | 3.3 |
| farthest 12, k = 2 | 7 | 300.0 | 1.4 | 10.2 | 2.5 | 2.4 |

## A4 random donor subsets (5 draws per target)

| donors | k | target × draw | Sn strict | Pr strict | F1 strict |
|:---|---:|---:|---:|---:|---:|
| 1 | 1 | 115 | 0.8 | 4.3 | 1.4 |
| 2 | 1 | 115 | 1.9 | 3.8 | 2.5 |
| 2 | 2 | 115 | 0.2 | 33.3 | 0.4 |
| 3 | 1 | 115 | 3.1 | 3.7 | 3.2 |
| 3 | 2 | 115 | 0.5 | 31.4 | 1.0 |
| 4 | 1 | 115 | 5.1 | 3.8 | 3.8 |
| 4 | 2 | 115 | 1.1 | 30.0 | 2.0 |
| 6 | 1 | 115 | 5.4 | 3.2 | 3.5 |
| 6 | 2 | 115 | 1.7 | 21.2 | 3.1 |
| 8 | 1 | 115 | 7.0 | 3.4 | 4.2 |
| 8 | 2 | 115 | 2.7 | 16.2 | 4.4 |
| 12 | 1 | 115 | 11.7 | 4.5 | 6.3 |
| 12 | 2 | 115 | 5.2 | 14.0 | 7.0 |
| 16 | 1 | 115 | 11.9 | 3.7 | 5.5 |
| 16 | 2 | 115 | 5.5 | 10.6 | 7.1 |

## A5 single-donor accuracy vs divergence (median de of the donor's projected mRNA)

| median de | (target, donor) pairs | lncRNA F1 (median) | mRNA F1 (median) |
|:---|---:|---:|---:|
| 0.00–0.02 | 14 | 27.7 | 50.2 |
| 0.02–0.04 | 28 | 5.9 | 33.9 |
| 0.04–0.06 | 40 | 2.7 | 24.4 |
| 0.06–0.10 | 98 | 2.7 | 24.6 |
| ≥ 0.10 | 326 | 0.7 | 11.1 |

## B pipeline configurations (lncRNA, fixed yardsticks, k = 2 unless stated)

| config | sequences projected | own reference (sum) | projected TD2 coding calls | raw models | of which single-exon | curated transcripts | Sn strict | Pr strict | F1 strict | ΔF1 vs B0 (median paired) | F1 lenient | F1 strict, k = 1 | chains = annotated mRNA | % chains on mRNA exon | curated decoys | decoy chains k = 1 | CPU-h re-executed |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| B0_baseline | 1,112,166 | 50,553 | 4,161 | 372,599 | 0 | 21,589 | 8.4 | 10.1 | 8.5 | +0.00 | 8.6 | 5.2 | 977 | 19.65 | 154 | 120,504 | 0.1 |
| B1_no_longest | 1,248,082 | 56,731 | 4,835 | 424,595 | 0 | 23,569 | 8.5 | 9.4 | 8.1 | -0.20 | 8.6 | 4.9 | 1,047 | 19.46 | 162 | 122,155 | 4.9 |
| B2_no_input_td2 | 1,129,788 | 51,354 | 10,326 | 375,905 | 0 | 21,713 | 8.5 | 10.0 | 8.5 | -0.01 | 8.6 | 5.2 | 997 | 19.85 | 154 | 120,504 | 2.5 |
| B3_no_proj_td2 | 1,112,166 | 50,553 | 0 | 376,760 | 0 | 21,722 | 8.4 | 10.0 | 8.4 | -0.03 | 8.6 | 5.2 | 1,050 | 19.99 | 154 | 120,504 | 0.3 |
| B4_no_td2 | 1,129,788 | 51,354 | 0 | 386,231 | 0 | 21,914 | 8.5 | 9.7 | 8.4 | -0.13 | 9.1 | 5.1 | 1,276 | 20.89 | 154 | 120,504 | 2.6 |
| B5_single_exon | 1,112,166 | 50,553 | 4,499 | 575,878 | 203,279 | 21,589 | 8.4 | 10.1 | 8.5 | +0.00 | 8.6 | 5.2 | 977 | 19.65 | 154 | 120,504 | 0.6 |
| B6_minimap2_ub | 1,112,166 | 50,553 | 4,186 | 379,679 | 0 | 21,683 | 8.4 | 10.1 | 8.4 | -0.05 | 8.6 | 5.1 | 978 | 19.63 | 170 | 125,147 | 7.9 |
| B7_minimap2_splicehq | 1,112,166 | 50,553 | 3,112 | 346,034 | 0 | 19,441 | 8.0 | 10.5 | 8.3 | +0.05 | 8.4 | 5.3 | 850 | 19.21 | 124 | 101,372 | 4.2 |

## Appendix — per target (lncRNA F1 strict)

| target | nearest donor (by median de) | its median de | all 22, k = 1 | all 22, k = 2 (default) | nearest 1, k = 1 | nearest 4, k = 2 | >= 2 congeners |
|:---|---:|---:|---:|---:|---:|---:|---:|
| Aricia agestis 91739 | Aricia artaxerxes 91738 | 0.007 | 9.2 | 26.8 | 46.1 | 39.9 | 1 |
| Aricia artaxerxes 91738 | Aricia agestis 91739 | 0.007 | 8.9 | 26.0 | 47.0 | 38.7 | 1 |
| Aricia cramera 282009 | Aricia artaxerxes 91738 | 0.009 | 7.6 | 25.4 | 37.6 | 37.4 | 1 |
| Aricia morronensis 1107572 | Aricia artaxerxes 91738 | 0.038 | 4.3 | 7.3 | 5.7 | 5.9 | 1 |
| Aricia nicias 689501 | Aricia artaxerxes 91738 | 0.021 | 2.3 | 6.2 | 4.2 | 4.9 | 1 |
| Celastrina argiolus 203782 | Glaucopsyche melanops 1002887 | 0.113 | 2.4 | 4.8 | 3.3 | 3.6 | 0 |
| Glaucopsyche melanops 1002887 | Glaucopsyche paphos 2744207 | 0.034 | 4.3 | 4.1 | 10.5 | 2.8 | 0 |
| Glaucopsyche paphos 2744207 | Glaucopsyche melanops 1002887 | 0.029 | 4.2 | 4.0 | 11.6 | 3.2 | 0 |
| Lycaena alciphron 282377 | Lycaena hippothoe 580924 | 0.036 | 13.6 | 30.4 | 29.2 | 37.1 | 1 |
| Lycaena hippothoe 580924 | Lycaena alciphron 282377 | 0.036 | 14.2 | 31.1 | 29.5 | 37.0 | 1 |
| Lycaena phlaeas 282391 | Lycaena alciphron 282377 | 0.051 | 7.7 | 12.8 | 11.4 | 14.2 | 1 |
| Lycaena thersamon 265360 | Lycaena virgaureae 282395 | 0.039 | 14.2 | 30.2 | 30.9 | 36.4 | 1 |
| Lycaena virgaureae 282395 | Lycaena alciphron 282377 | 0.038 | 14.4 | 32.5 | 31.1 | 38.9 | 1 |
| Lysandra albicans 282008 | Lysandra coridon 268709 | 0.008 | 7.5 | 17.9 | 25.5 | 23.7 | 1 |
| Lysandra bellargus 138070 | Lysandra albicans 282008 | 0.009 | 8.2 | 17.5 | 26.2 | 23.4 | 1 |
| Lysandra coridon 268709 | Lysandra albicans 282008 | 0.008 | 7.7 | 17.9 | 26.0 | 24.2 | 1 |
| Plebejus argus 242267 | Plebejus idas 155166 | 0.047 | 2.6 | 5.2 | 4.0 | 3.7 | 0 |
| Plebejus idas 155166 | Plebejus argus 242267 | 0.054 | 2.5 | 4.5 | 2.8 | 2.6 | 0 |
| Polyommatus atlantica 712039 | Polyommatus icarus 265386 | 0.035 | 1.8 | 5.2 | 3.1 | 3.8 | 1 |
| Polyommatus icarus 265386 | Polyommatus szabokyi 2079168 | 0.008 | 5.2 | 8.5 | 11.7 | 4.5 | 1 |
| Polyommatus szabokyi 2079168 | Polyommatus icarus 265386 | 0.008 | 5.2 | 3.7 | 9.2 | 2.0 | 1 |
| Satyrium ilicis 876074 | Satyrium pruni 876075 | 0.053 | 3.2 | 2.1 | 5.7 | 1.4 | 0 |
| Satyrium pruni 876075 | Satyrium ilicis 876074 | 0.060 | 4.8 | 2.9 | 6.7 | 1.5 | 0 |
