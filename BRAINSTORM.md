This file contains the high level descripion of the fomo pipleine. 

The goal of the pipeline is to take a target .fa assembly as input and return a .gff3 of candidate long non-coding RNA annotated on the target assembly.
The candidates long non-coding RNA are identified in the target assembly using lncRNA annotation from source assembly/annotations of other species.  

The high level steps are:
- collect source annotation (defer implementation)
- preprocessing
- projection
- validation
- benchmarking (only if reference annotation is available)
- reporting

### Next step to implememt
- add samtools stats after the alignement step and send the result to MultiQC
