#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { get_bed_ref }                from '../modules/utils.nf'
include { publish as publishConsensus} from '../modules/utils.nf'
include { publish as publishQCCSV}     from '../modules/utils.nf'
include { performHostFilter }          from '../modules/illumina.nf'
include { align_trim }                 from '../modules/illumina.nf'
include { indexReference}              from '../modules/illumina.nf'
include { readMapping }                from '../modules/illumina.nf'
include { callConsensusFreebayes }     from '../modules/illumina.nf'
include { annotateVariantsVCF }        from '../modules/illumina.nf'
include { alignConsensusToReference }  from '../modules/illumina.nf'

include { makeQCCSV }         from '../modules/qc.nf'
// include { writeQCSummaryCSV } from '../modules/qc.nf'

include { collateSamples }    from '../modules/upload.nf'

workflow prepareReferenceFiles {

  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";


  // Scheme selection logic grabbed from epi2melabs/wf-artic, thanks ONT!
  if (!params.custom_scheme){

    schemes = file(projectDir.resolve("./data/primer-schemes/**bed"), type: 'file', maxdepth: 10)

    valid_scheme_versions = []

    log.info """
    ------------------------------------
    Available Primer Schemes:
    ------------------------------------
    """
    log.info """  Name\t\tVersion"""
    for (scheme in schemes){
      main = scheme.toString().split("primer-schemes/")[1]
      name = main.split("/")[0]
      version = """${main.split("/")[1]}/${main.split("/")[2]}"""
      valid_scheme_versions.add(version)
      log.info """${c_green}  ${name}\t${version}\t${c_reset}"""
    }

    log.info """
    ------------------------------------
    """

    if (params.list_schemes) {
      exit 1
    }

    if (!valid_scheme_versions.any { it == params.scheme_version}) {
        println("`--scheme_version` should be one of: $valid_scheme_versions, for `--scheme_name`: $params.scheme_name")
        exit 1
    }
 
    scheme_dir_name = "primer-schemes"
    schemes = """./data/${scheme_dir_name}/${params.scheme_name}"""
    scheme_dir = file(projectDir.resolve(schemes), type:'file', checkIfExists:true)

    primers_path = """./data/${scheme_dir_name}/${params.scheme_name}/${params.scheme_version}/primer.bed"""
    primers = file(projectDir.resolve(primers_path), type:'file', checkIfExists:true)

    reference_path = """./data/${scheme_dir_name}/${params.scheme_name}/${params.scheme_version}/reference.fasta"""
    reference = file(projectDir.resolve(reference_path),type:'file', checkIfExists:true)

    params._scheme_version = params.scheme_version
    params._scheme_name = params.scheme_name

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
        performHostFilter(ch_filePairs)
        ch_filtered_reads = performHostFilter.out
      } else {
        ch_filtered_reads = ch_filePairs
      }

      readMapping(ch_filtered_reads.combine(ch_preparedRef), ch_bwaIndex)

      align_trim(readMapping.out, ch_bedFile)

      callConsensusFreebayes(align_trim.out.ptrimmed_bam.combine(ch_preparedRef))

      if (params.align_consensus) {
        alignConsensusToReference(callConsensusFreebayes.out.consensus.combine(ch_preparedRef))
        alignConsensusToReference.out.map{ sampleName,sampleFasta -> sampleFasta }.collectFile(name: "all_consensus.aln.fa").set{ alignment }

        publishConsensus(alignment)
      }
      

      makeQCCSV(align_trim.out.ptrimmed_bam.join(callConsensusFreebayes.out.consensus, by: 0)
          .combine(ch_preparedRef)
				  .combine(ch_bedFile)
          )

      makeQCCSV.out.csv
          .collectFile(name: "${params.prefix}.qc.csv", skip: 1, keepHeader: true).set { qc }
      
      publishQCCSV(qc)

    emit:
      qc_pass = callConsensusFreebayes.out.consensus.join(ch_filtered_reads)
}

workflow mpxvIllumina {
    take:
      ch_filePairs

    main:
      prepareReferenceFiles()
      sequenceAnalysis(ch_filePairs, prepareReferenceFiles.out.reference, prepareReferenceFiles.out.bwaIndex, prepareReferenceFiles.out.primers)
}


