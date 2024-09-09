#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { get_bed_ref }                from '../modules/utils.nf'
include { publish as publishConsensus} from '../modules/utils.nf'
include { publish as publishQCCSV}     from '../modules/utils.nf'
include { normalizeDepth }             from '../modules/illumina.nf'
include { performHostFilter }          from '../modules/illumina.nf'
include { readTrimming }               from '../modules/illumina.nf'
include { filterResidualAdapters }     from '../modules/illumina.nf'
include { indexReference}              from '../modules/illumina.nf'
include { readMapping }                from '../modules/illumina.nf'
include { trimPrimerSequences }        from '../modules/illumina.nf'
include { callConsensusFreebayes }     from '../modules/illumina.nf'
include { annotateVariantsVCF }        from '../modules/illumina.nf'
include { alignConsensusToReference }  from '../modules/illumina.nf'

include { makeQCCSV }         from '../modules/qc.nf'
// include { writeQCSummaryCSV } from '../modules/qc.nf'

include { collateSamples }    from '../modules/upload.nf'

include { squirrelAlignmentAndQC } from '../modules/squirrel.nf'

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

      if (!params.skip_normalize_depth) {
        ch_reads_to_hostfilter = normalizeDepth(ch_filePairs)
      } else {
        ch_reads_to_hostfilter = ch_filePairs
      }

      performHostFilter(ch_reads_to_hostfilter)

      readTrimming(performHostFilter.out.fastqPairs)

      filterResidualAdapters(readTrimming.out)

      readMapping(filterResidualAdapters.out.combine(ch_preparedRef), ch_bwaIndex)

      trimPrimerSequences(readMapping.out.combine(ch_bedFile))

      callConsensusFreebayes(trimPrimerSequences.out.ptrim.combine(ch_preparedRef))

      if (params.gff != 'NO_FILE') {
        annotateVariantsVCF(callConsensusFreebayes.out.variants.combine(ch_preparedRef).combine(ch_gff))
      }

      alignConsensusToReference(callConsensusFreebayes.out.consensus.combine(ch_preparedRef))
      alignConsensusToReference.out.map{ sampleName,sampleFasta -> sampleFasta }.collectFile(name: "all_consensus.aln.fa").set{ alignment }

      publishConsensus(alignment)

      makeQCCSV(trimPrimerSequences.out.ptrim.join(callConsensusFreebayes.out.consensus, by: 0)
          .combine(ch_preparedRef)
				  .combine(ch_bedFile)
          )

      makeQCCSV.out.csv
          .collectFile(name: "${params.prefix}.qc.csv", skip: 1, keepHeader: true).set { qc }
      
      publishQCCSV(qc)

      collateSamples(callConsensusFreebayes.out.consensus.join(performHostFilter.out.fastqPairs))

      callConsensusFreebayes.out.consensus.map{ sampleName,sampleFasta -> sampleFasta }.collectFile(name: "all_consensus.fa").set{ consensus }
      if ( params.squirrel_assembly_refs ) {
            refs_ch = channel.fromPath("${params.squirrel_assembly_refs}", checkIfExists:true)
      } else {
            refs_ch = channel.fromPath("${projectDir}/test_data/empty.fa", checkIfExists:true)
      }
      squirrelAlignmentAndQC(consensus, refs_ch)

    emit:
      qc_pass = collateSamples.out
      alignment = squirrelAlignmentAndQC.out.alignment
}

workflow mpxvIllumina {
    take:
      ch_filePairs

    main:
      prepareReferenceFiles()
      sequenceAnalysis(ch_filePairs, prepareReferenceFiles.out.reference, prepareReferenceFiles.out.bwaIndex, prepareReferenceFiles.out.bedfile, prepareReferenceFiles.out.gff)


}


