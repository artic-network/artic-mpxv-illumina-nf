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
    // Get reference fasta
    // if (params.ref) {
    //   Channel.fromPath(params.ref)
    //           .set{ ch_refFasta }
    // } else {
    //   articDownloadScheme()
    //   articDownloadScheme.out.reffasta
    //                       .set{ ch_refFasta }
    // }

    scheme_dir_name = "primer-schemes"
    schemes = """./data/${scheme_dir_name}/${params.scheme_name}"""
    scheme_dir = file(projectDir.resolve(schemes), type:'file', checkIfExists:true)

    get_bed_ref(scheme_dir, params.scheme_name, params.scheme_version)
    ch_bedFile = get_bed_ref.out.bed
    ch_refFasta = get_bed_ref.out.ref
    ch_gff = get_bed_ref.out.gff


    /* Either get BWA aux files from reference 
      location or make them fresh */
  
    // Index the reference
    indexReference(ch_refFasta)

    indexReference.out
      .set{ ch_bwaIndex }

    emit:
      bwaIndex = ch_bwaIndex
      reference = ch_refFasta
      bedfile = ch_bedFile
      gff = ch_gff
}


workflow sequenceAnalysis {
    take:
      ch_filePairs
      ch_preparedRef
      ch_bwaIndex
      ch_bedFile
      ch_gff

    main:

      // if (!params.skip_normalize_depth) {
      //   ch_reads_to_hostfilter = normalizeDepth(ch_filePairs)
      // } else {
      //   ch_reads_to_hostfilter = ch_filePairs
      // }

      performHostFilter(ch_filePairs)

      readMapping(performHostFilter.out.combine(ch_preparedRef), ch_bwaIndex)

      align_trim(readMapping.out, ch_bedFile)

      callConsensusFreebayes(align_trim.out.ptrimmed_bam.combine(ch_preparedRef))

      if (params.gff != 'NO_FILE') {
        annotateVariantsVCF(callConsensusFreebayes.out.variants.combine(ch_preparedRef).combine(ch_gff))
      }

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
      qc_pass = callConsensusFreebayes.out.consensus.join(performHostFilter.out)
}

workflow mpxvIllumina {
    take:
      ch_filePairs

    main:
      prepareReferenceFiles()
      sequenceAnalysis(ch_filePairs, prepareReferenceFiles.out.reference, prepareReferenceFiles.out.bwaIndex, prepareReferenceFiles.out.bedfile, prepareReferenceFiles.out.gff)
}


