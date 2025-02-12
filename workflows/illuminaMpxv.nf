#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fetchHostileReference }      from '../modules/illumina.nf'
include { performHostFilter }          from '../modules/illumina.nf'
include { align_trim }                 from '../modules/illumina.nf'
include { readTrimming }               from '../modules/illumina.nf'
include { readMapping }                from '../modules/illumina.nf'
include { callConsensusFreebayes }     from '../modules/illumina.nf'
include { annotateVariantsVCF }        from '../modules/illumina.nf'
include { alignConsensusToReference }  from '../modules/illumina.nf'
include { fetchScheme }                from '../modules/illumina.nf'
include { fetchScheme as fetchScheme_free } from '../modules/illumina.nf'

include { makeQCCSV }         from '../modules/qc.nf'

include { squirrelAlignmentAndQC } from '../modules/squirrel.nf'

workflow prepareReferenceFiles {

  take:
    ch_trimmed_reads
  
  main:

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";


    // Scheme selection logic grabbed from epi2melabs/wf-artic, thanks ONT!
    if (!params.custom_scheme){

      if (!params.freetext_scheme_name) {
        fetchScheme(params.scheme_name, ch_trimmed_reads)

        ch_scheme = fetchScheme.out
      } else {
        fetchScheme_free(params.freetext_scheme_name, ch_trimmed_reads)

        ch_scheme = fetchScheme_free.out
      }

    } else {
      //custom scheme path defined
      log.info """${c_purple}Custom primer scheme selected: ${params.custom_scheme} (WARNING: We do not validate your scheme - use at your own risk!)${c_reset}"""
      //check path for re§ired files
      primers = file("""${params.custom_scheme}/primer.bed""", type:'file', checkIfExists:true)
      reference = file("""${params.custom_scheme}/reference.fasta""", type:'file', checkIfExists:true)

      scheme = Channel.of([reference, primers])

      ch_scheme = ch_trimmed_reads
        .map { sample, forward, reverse -> sample}
        .combine(scheme)

    }    

    emit:
      ch_scheme
}


workflow sequenceAnalysis {
    take:
      ch_filePairs

    main:

      if (!params.skip_host_filter) {
        fetchHostileReference()
        performHostFilter(ch_filePairs, fetchHostileReference.out)
        ch_filtered_reads = performHostFilter.out
      } else {
        ch_filtered_reads = ch_filePairs
      }

      readTrimming(ch_filtered_reads)

      prepareReferenceFiles(ch_filtered_reads)

      readTrimming.out.combine(prepareReferenceFiles.out, by: 0)
        .set{ ch_toMap }

      readMapping(ch_toMap)

      align_trim(readMapping.out)

      callConsensusFreebayes(align_trim.out.ptrimmed_bam)

      if (params.align_consensus) {
        alignConsensusToReference(callConsensusFreebayes.out.consensus)
        alignConsensusToReference.out
          .map{ sampleName, sampleFasta -> sampleFasta }
          .collectFile(name: "${params.prefix}.all_consensus.aln.fasta", storeDir: "${params.outdir}/")
          .set{ alignment }
      }

      ch_qc = align_trim.out.ptrimmed_bam
        .combine(callConsensusFreebayes.out.consensus, by: 0)
      
      makeQCCSV(ch_qc)

      makeQCCSV.out.csv
          .collectFile(name: "${params.prefix}.qc.csv", skip: 1, keepHeader: true, storeDir: "${params.outdir}/")
          .set { qc }
      
      callConsensusFreebayes.out.consensus.map{ sampleName, sampleFasta -> sampleFasta }
        .collectFile(name: "${params.prefix}.all_consensus.fasta", storeDir: "${params.outdir}/")
        .set{ consensus }

      if ( params.squirrel_assembly_refs ) {
            refs_ch = channel.fromPath("${params.squirrel_assembly_refs}", checkIfExists:true)
      } else {
            refs_ch = channel.fromPath("${projectDir}/test_data/empty.fasta", checkIfExists:true)
      }
      
      if (!params.skip_squirrel) {
        squirrelAlignmentAndQC(consensus, refs_ch)
      }

}

