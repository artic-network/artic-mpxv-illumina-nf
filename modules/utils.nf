process publish {
    publishDir "${params.outdir}/", mode: 'copy'
    container 'nextflow/bash:latest'

    input:
        path name
    output:
        path name
    script:
    """
    """
}