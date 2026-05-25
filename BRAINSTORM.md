This file contains the high level descripion of the fomo pipleine. 

The goal of the pipeline is to take a target .fa assembly as input and return a .gff3 of candidate long non-coding RNA annotated on the target assembly.
The candidates long non-coding RNA are identified in the target assembly using lncRNA annotation from source assembly/annotations of other species.  

The high level steps are:
- collect source annotation
- preprocessing
- projection
- validation
- reporting


### preprocessing
- extract lncRNA and mRNA spliced fasta sequence (legacy_script/preprocessing/00_get_feature_annotation.sh)
- extract intergenic intervalas (legacy_script/preprocessing/02_get_intergenic_intervals.sh)
- build decoy sequence (legacy_script/preprocessing/03_relocate_loci.py)
- extract decoy spliced fasta sequence (legacy_script/preprocessing/04_get_decoy_sequence.sh
- build preporcessing statistics (legacy_script/preprocessing/05_get_gff_statistics.sh)

### projection (mapping)
- map fasta sequence (legacy_script/minimap_transfer/00_run_minimap_base.sh)
- convert bam to gff (legacy_script/minimap_transfer/01_convert_bam_to_gff.sh)


### statistics 
- preprocessing statistics
- extract alignment metrics
- extract projected gff metrics
