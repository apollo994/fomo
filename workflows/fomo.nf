include { PREPROCESSING } from '../subworkflows/local/preprocessing'

workflow FOMO {
    // Parse samplesheet
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            assert row.species : "Samplesheet row is missing 'species': ${row}"
            assert row.role in ['source', 'target'] : "Invalid role '${row.role}' for ${row.species} (must be 'source' or 'target')"
            assert row.fasta : "Samplesheet row for ${row.species} is missing 'fasta'"
            assert row.gff3  : "Samplesheet row for ${row.species} is missing 'gff3'"

            def meta  = [id: row.species, role: row.role]
            def fasta = row.fasta.startsWith('/') ? file(row.fasta, checkIfExists: true)
                                                   : file("${projectDir}/${row.fasta}", checkIfExists: true)
            def gff3  = row.gff3.startsWith('/') ? file(row.gff3, checkIfExists: true)
                                                  : file("${projectDir}/${row.gff3}", checkIfExists: true)
            tuple(meta, fasta, gff3)
        }
        .branch {
            target: it[0].role == 'target'
            source: it[0].role == 'source'
        }
        .set { ch_input }

    PREPROCESSING(ch_input.source)
}
