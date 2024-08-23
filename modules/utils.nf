process articDownloadScheme{
    tag params.schemeRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scheme", mode: "copy"

    output:
    path "scheme/**/${params.schemeVersion}/*.reference.fasta" , emit: reffasta
    path "scheme/**/${params.schemeVersion}/${params.scheme}.bed" , emit: bed
    path "scheme" , emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} scheme
    """
}

process get_bed_ref {

    label 'process_single'

    container 'nextflow/bash:latest'

    input:
        path scheme_dir
        val scheme_name
        val scheme_version
    output:
        path "primer.bed", emit: bed
        path "reference.fasta", emit: ref
        path "reference.gff3", emit: gff

    """
    cp ${scheme_name}/${scheme_version}/primer.bed primer.bed
    cp ${scheme_name}/${scheme_version}/reference.fasta reference.fasta
    cp ${scheme_name}/${scheme_version}/reference.gff3 reference.gff3
    """
}

process performHostFilter {

    tag { sampleName }

    container 'community.wave.seqera.io/library/bwa_pysam_samtools_python:f3b7e3fe2ad2cadc'

    conda 'bioconda::bwa=0.7.18', 'bioconda::samtools=1.20', 'bioconda::python=3.12.5', 'bioconda::pysam=0.22.1'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_hostfiltered_R*.fastq.gz", mode: 'copy'

    input:
        tuple val(sampleName), path(forward), path(reverse)
    output:
        tuple val(sampleName), path("${sampleName}_hostfiltered_R1.fastq.gz"), path("${sampleName}_hostfiltered_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        bwa mem -t ${task.cpus} ${params.composite_ref} ${forward} ${reverse} | \
            filter_non_human_reads.py -c ${params.viral_contig_name} > ${sampleName}.viral_and_nonmapping_reads.bam
        samtools sort -n ${sampleName}.viral_and_nonmapping_reads.bam | \
             samtools fastq -1 ${sampleName}_hostfiltered_R1.fastq.gz -2 ${sampleName}_hostfiltered_R2.fastq.gz -s ${sampleName}_singletons.fastq.gz -
        """
}

process publish {
    publishDir "${params.outdir}/", mode: 'copy'
    input:
        path name
    output:
        path name
    script:
    """
    """
}