#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fetchHostileReference }      from '../modules/illumina.nf'
include { performHostFilter }          from '../modules/illumina.nf'
include { align_trim }                 from '../modules/illumina.nf'
include { indexReference}              from '../modules/illumina.nf'
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

  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";


  // Scheme selection logic grabbed from epi2melabs/wf-artic, thanks ONT!
  if (!params.custom_scheme){

    if (!params.freetext_scheme_name) {
      fetchScheme(params.scheme_name)

      reference = fetchScheme.out.reference
      primers = fetchScheme.out.primers
    } else {
      fetchScheme_free(params.freetext_scheme_name)

      reference = fetchScheme_free.out.reference
      primers = fetchScheme_free.out.primers
    }

  } else {
    //custom scheme path defined
    log.info """${c_purple}Custom primer scheme selected: ${params.custom_scheme} (WARNING: We do not validate your scheme - use at your own risk!)${c_reset}"""
    //check path for re§ired files
    primers = file("""${params.custom_scheme}/primer.bed""", type:'file', checkIfExists:true)
    reference = file("""${params.custom_scheme}/reference.fasta""", type:'file', checkIfExists:true)

    params._scheme_version = 'None'
    params._scheme_name = params.scheme_name

    scheme_dir =  params.custom_scheme
  }    

    scheme_dir_name = "primer-schemes"
    schemes = """./data/${scheme_dir_name}/${params.scheme_name}"""
    scheme_dir = file(projectDir.resolve(schemes), type:'file', checkIfExists:true)

    /* Either get BWA aux files from reference 
      location or make them fresh */
  
    // Index the reference
    indexReference(reference)

    indexReference.out
      .set{ ch_bwaIndex }

    emit:
      bwaIndex = ch_bwaIndex
      reference
      primers
}


workflow sequenceAnalysis {
    take:
      ch_filePairs
      ch_preparedRef
      ch_bwaIndex
      ch_bedFile

    main:

      if (!params.skip_host_filter) {
        fetchHostileReference()
        performHostFilter(ch_filePairs, fetchHostileReference.out)
        ch_filtered_reads = performHostFilter.out
      } else {
        ch_filtered_reads = ch_filePairs
      }

      readTrimming(ch_filtered_reads)

      readMapping(readTrimming.out.combine(ch_preparedRef), ch_bwaIndex)

      align_trim(readMapping.out, ch_bedFile)

      callConsensusFreebayes(align_trim.out.ptrimmed_bam.combine(ch_preparedRef))

      if (params.align_consensus) {
        alignConsensusToReference(callConsensusFreebayes.out.consensus.combine(ch_preparedRef))
        alignConsensusToReference.out
          .map{ sampleName,sampleFasta -> sampleFasta }
          .collectFile(name: "${params.prefix}.all_consensus.aln.fasta", storeDir: "${params.outdir}/")
          .set{ alignment }
      }
      
      makeQCCSV(align_trim.out.ptrimmed_bam.join(callConsensusFreebayes.out.consensus, by: 0)
          .combine(ch_preparedRef)
				  .combine(ch_bedFile)
          )

      makeQCCSV.out.csv
          .collectFile(name: "${params.prefix}.qc.csv", skip: 1, keepHeader: true, storeDir: "${params.outdir}/")
          .set { qc }
      
      callConsensusFreebayes.out.consensus.map{ sampleName,sampleFasta -> sampleFasta }
        .collectFile(name: "${params.prefix}.all_consensus.fasta", storeDir: "${params.outdir}/")
        .set{ consensus }

      if ( params.squirrel_assembly_refs ) {
            refs_ch = channel.fromPath("${params.squirrel_assembly_refs}", checkIfExists:true)
      } else {
            refs_ch = channel.fromPath("${projectDir}/test_data/empty.fasta", checkIfExists:true)
      }
      squirrelAlignmentAndQC(consensus, refs_ch)

    emit:
      alignment = squirrelAlignmentAndQC.out.alignment
}

workflow mpxvIllumina {
    take:
      ch_filePairs

    main:
      prepareReferenceFiles()

      sequenceAnalysis(ch_filePairs, prepareReferenceFiles.out.reference, prepareReferenceFiles.out.bwaIndex, prepareReferenceFiles.out.primers)
}


