#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FOMO } from './workflows/fomo'

workflow {
    FOMO()
}
