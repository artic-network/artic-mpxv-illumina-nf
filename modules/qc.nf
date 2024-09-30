process makeQCCSV {

    tag { sampleName }

    label 'process_single'

    container 'community.wave.seqera.io/library/samtools_biopython_matplotlib_pandas:2bb143fe29dc1ce9'

    conda "conda-forge::matplotlib=3.9.2 conda-forge::pandas=2.2.2 conda-forge::biopython=1.84 bioconda::samtools=1.20"

    publishDir "${params.outdir}/qc_plots", pattern: "${sampleName}.depth.png", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bam_index), path(fasta), path(ref), path(primer_bed)

    output:
    path "${params.prefix}.${sampleName}.qc.csv", emit: csv
    path "${sampleName}.depth.png"

    script:
    """
    qc.py --outfile ${params.prefix}.${sampleName}.qc.csv --sample ${sampleName} --ref ${ref} --bam ${bam} --fasta ${fasta} --primer-bed ${primer_bed} --min-depth ${params.varMinDepth} --primer-pairs primer-pairs.tsv
    """
}


process writeQCSummaryCSV {
    tag { params.prefix }

    label 'process_single'

    container 'jitesoft/alpine:3.20.2'

    input:
    val lines

    exec:
    file("${params.outdir}/${params.prefix}.qc.csv").withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
}

