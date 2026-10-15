# Plans

Design decisions, one file per plan. Plans 00–26 carry a decision-record header (status, date,
decision); their bodies are kept as written. Plan 27 is the refactor that produced the current code
base: a master plan with subplans 27a–27k, a frozen copy of the plan as first approved in
[27_refactor_original/](27_refactor_original/), and the review of every deviation in
`27z_deviation_review.md` once 27k is done.

| plan | title | status | date |
|:---|:---|:---|:---|
| [00_preprocessing-step00.md](00_preprocessing-step00.md) | Plan: FOMO Preprocessing — Step 00 (incremental) | Partially superseded by 20 and 27 | 2026-05-24 |
| [01_preprocessing-decoy.md](01_preprocessing-decoy.md) | Plan: FOMO Preprocessing — Decoy track (steps 02–04) | Implemented (optional since 12) | 2026-05-24 |
| [02_input-validation.md](02_input-validation.md) | Input validation with nf-schema | Implemented | 2026-05-24 |
| [03_mapping.md](03_mapping.md) | Plan: FOMO Mapping — Project spliced FASTAs onto target | Superseded by 15 | 2026-05-25 |
| [04_bam_to_gff.md](04_bam_to_gff.md) | Plan: FOMO BAM → GFF3 conversion | Implemented | 2026-05-25 |
| [05_unify_projection.md](05_unify_projection.md) | Plan: Unify MAPPING + BAM_TO_GFF under PROJECTION + record alignment % | Implemented; partially superseded by 15 | 2026-05-25 |
| [06_benchmarking.md](06_benchmarking.md) | Plan: FOMO Benchmarking — gffcompare per projection | Partially superseded by 18 and 20 | 2026-05-26 |
| [07_reporting_refactor.md](07_reporting_refactor.md) | Plan: refactor STATS → REPORTING, stage-local QC, split MultiQC config | Partially superseded by 16 and 27 | 2026-05-26 |
| [08_mixed_compression_inputs.md](08_mixed_compression_inputs.md) | Plan: Accept both compressed (.gz) and uncompressed FASTA / GFF3 inputs | Superseded by 26 | 2026-05-27 |
| [09_top3_consensus.md](09_top3_consensus.md) | Plan: Top-3 consensus annotation (selected by lncRNA transcript F1) | Superseded by 23 | 2026-05-27 |
| [10_replace_agat_stats.md](10_replace_agat_stats.md) | Plan: replace AGAT statistics with bin/gff_stats_standalone.py | Superseded by 11 | 2026-09-22 |
| [11_rust_gff_feature_stats.md](11_rust_gff_feature_stats.md) | Plan: replace AGAT statistics with the Rust gff-feature-stats | Implemented | 2026-09-22 |
| [12_optional_mrna_decoy.md](12_optional_mrna_decoy.md) | Plan: make the mRNA and decoy tracks optional (--include_mrna / --include_decoy) | Implemented | 2026-09-22 |
| [13_resource_rightsizing.md](13_resource_rightsizing.md) | Plan: Resource right-sizing from execution-trace evidence | Partially superseded by 24 and 27 | 2026-09-22 |
| [14_multitarget.md](14_multitarget.md) | 14 — Multi-target: role: both and all-vs-all in one run | Implemented | 2026-09-22 |
| [15_single_mapping.md](15_single_mapping.md) | 15 — Single mapping: one all-sources alignment per target | Implemented | 2026-09-22 |
| [15a_single_mapping_align.md](15a_single_mapping_align.md) | 15a — Phase A: merge, align once, split back | Implemented | 2026-09-22 |
| [15b_single_mapping_aggregates.md](15b_single_mapping_aggregates.md) | 15b — Phase B: raw and collapsed aggregates | Partially superseded by 23 | 2026-09-22 |
| [15c_single_mapping_outputs.md](15c_single_mapping_outputs.md) | 15c — Phase C: output layout, reporting, resources, docs | Implemented | 2026-09-22 |
| [16_general_report.md](16_general_report.md) | 16 — Run-level "general" MultiQC report | Implemented | 2026-09-22 |
| [17_scheduling_efficiency.md](17_scheduling_efficiency.md) | Plan: Scheduling efficiency — executor tuning + removing redundant gunzip tasks | Implemented; partially superseded by 26 | 2026-09-22 |
| [18_per_target_batching.md](18_per_target_batching.md) | 18 — Per-target batching: one process for all sources | Implemented; partially superseded by 27 | 2026-09-08 |
| [19_run_summary_scaling.md](19_run_summary_scaling.md) | 19 — RUN_SUMMARY_TABLES at scale: why it timed out, and what to do about it | Superseded by 27 | 2026-09-22 |
| [20_longest_isoform.md](20_longest_isoform.md) | 20 — One filtered annotation per species: longest isoform, used for projection, reference and results | Implemented | 2026-10-01 |
| [21_curated_annotation.md](21_curated_annotation.md) | 21 — Curated annotation: chains shared by ≥ 2 species, merged with the reference | Implemented | 2026-10-06 |
| [22_annotated_targets.md](22_annotated_targets.md) | 22 — Annotated targets: a target may carry a gff3, used as the curation reference | Implemented | 2026-10-07 |
| [23_no_self_no_top3.md](23_no_self_no_top3.md) | 23 — Remove self-projection, the top-N consensus and allModels_collapsed | Implemented | 2026-10-07 |
| [24_resource_rightsizing.md](24_resource_rightsizing.md) | Plan: CPU/memory right-sizing — no CPU escalation on retry, smaller batch and minimap2 requests | Implemented; partially superseded by 27 | 2026-10-08 |
| [25_curated_model_format.md](25_curated_model_format.md) | 25 — Curated GFF3 format: content-hashed IDs, one transcript per supported chain | Implemented | 2026-10-08 |
| [26_compressed_genome_io.md](26_compressed_genome_io.md) | Plan: keep genomes compressed — decompress on node-local disk, never on NFS | Implemented | 2026-10-14 |
| [27_refactor.md](27_refactor.md) | 27 — Codebase refactor: reduce complexity (master plan) | in progress (started 2026-10-15) |  |
| [27a_safety_net.md](27a_safety_net.md) | 27a — Safety net: golden baseline, pytest, nf-test, lint, CI | done (2026-10-15) |  |
| [27b_delete_and_fix.md](27b_delete_and_fix.md) | 27b — Delete dead code, fix latent bugs, comment sweep | done (2026-10-15) |  |
| [27c_local_modules.md](27c_local_modules.md) | 27c — All modules local; containers pinned in one place; versions collected; no conda | done (2026-10-15) |  |
| [27d_config_structure.md](27d_config_structure.md) | 27d — Config structure: include order, selectors, resources, site profile, pipeline_info | done (2026-10-15) |  |
| [27e_track_identity.md](27e_track_identity.md) | 27e — Track registry, meta identity, channel simplification | done (2026-10-15) |  |
| [27f_python_shared.md](27f_python_shared.md) | 27f — Shared Python modules, filter_transcript.py, curate_models split | done (2026-10-15) |  |
| [27g_stats_records.md](27g_stats_records.md) | 27g — Producers emit data records; per-source identity from a manifest | done (2026-10-15) |  |
| [27h_run_summary_records.md](27h_run_summary_records.md) | 27h — The run summary reads data records | done (2026-10-15) |  |
| [27i_multiqc_from_records.md](27i_multiqc_from_records.md) | 27i — MultiQC rendered from records; adapters, duplicates and the old TSV path deleted | done (2026-10-15) |  |
| [27j_docs.md](27j_docs.md) | 27j — Documentation, plan records, dev tools, comment budget | done (2026-10-15) |  |
| [27k_final_validation.md](27k_final_validation.md) | 27k — Final validation, real-data regression, plan comparison | not started |  |
