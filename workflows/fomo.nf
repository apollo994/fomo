include { PREPROCESSING     } from '../subworkflows/local/preprocessing'
include { STATS             } from '../subworkflows/local/stats'
include { PROJECTION        } from '../subworkflows/local/projection'
include { samplesheetToList } from 'plugin/nf-schema'

workflow FOMO {
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .branch { meta, fasta, gff3 ->
            target: meta.role == 'target'
            source: meta.role == 'source'
        }
        .set { ch_input }

    PREPROCESSING(ch_input.source)

    STATS(
        ch_input.source,
        ch_input.target,
        PREPROCESSING.out.filtered_gff3,
        PREPROCESSING.out.decoy_gff3,
        PREPROCESSING.out.spliced_fasta,
        PREPROCESSING.out.decoy_spliced_fasta
    )

    PROJECTION(
        ch_input.target,
        PREPROCESSING.out.spliced_fasta,
        PREPROCESSING.out.decoy_spliced_fasta
    )
}
