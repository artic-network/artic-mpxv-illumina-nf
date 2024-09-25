process squirrelAlignmentAndQC {

    label "process_low"

    container "docker.io/articnetworkorg/squirrel@sha256:78921da0c726dd525f5d72db331526975adb634478ee31047c46765bd3d5a64a"

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "squirrel/**", mode: 'copy', saveAs: {  fn -> fn.replace("squirrel/", "")}

    input:
        path fasta
        path refs
    output:
        path "squirrel/all_consensus.aln.fasta", emit: alignment
        path "squirrel/**", emit: all
        path "squirrel.version", emit: version

  script:
    extra = ""
    if ( params.squirrel_assembly_refs )
        extra += " --assembly-refs ${refs}"
    if ( params.clade )
        extra += " --clade ${params.clade}"
    if ( params.run_phylo )
        extra += " --run-phylo"
    if ( params.run_seq_qc )
        extra += " --run-seq-qc"
    if ( params.outgroups )
        extra += " --outgroups ${params.outgroups}"

    """
    export XDG_CACHE_HOME=\$PWD/.cache
    squirrel --version 2>&1 | sed 's/: /,/' > squirrel.version
    squirrel ${fasta} --no-mask --seq-qc -o squirrel --outfile all_consensus.aln.fasta --tempdir squirrel_tmp -t ${task.cpus} ${extra}
    """
}
