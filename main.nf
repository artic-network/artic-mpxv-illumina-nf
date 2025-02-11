#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import subworkflows
include {sequenceAnalysis} from './workflows/illuminaMpxv.nf'

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)

// main workflow
workflow {

  Channel.fromFilePairs( params.fastqSearchPath, flat: true)
          .filter{ !( it[0] =~ /Undetermined/ ) }
          .set{ ch_filePairs }

  main:
    sequenceAnalysis(ch_filePairs)
     
}
