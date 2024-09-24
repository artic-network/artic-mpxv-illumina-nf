#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// include modules
include {printHelp} from './modules/help.nf'

// import subworkflows
include {mpxvIllumina} from './workflows/illuminaMpxv.nf'

if (params.help){
    printHelp()
    exit 0
}

if ( !params.directory ) {
  println("Please supply a directory containing fastqs or CRAMs with --directory.")
  println("Use --help to print help")
  System.exit(1)
}

if ( ! params.prefix ) {
     println("Please supply a prefix for your output files with --prefix")
     println("Use --help to print help")
     System.exit(1)
} else {
     if ( params.prefix =~ /\// ){
         println("The --prefix that you supplied contains a \"/\", please replace it with another character")
         System.exit(1)
     }
} 

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)

// main workflow
workflow {

  Channel.fromFilePairs( params.fastqSearchPath, flat: true)
          .filter{ !( it[0] =~ /Undetermined/ ) }
          .set{ ch_filePairs }

  main:
    mpxvIllumina(ch_filePairs)
     
}
