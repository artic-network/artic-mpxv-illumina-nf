process squirrelAlignmentAndQC {

    label "process_low"

    container "docker.io/articnetworkorg/squirrel@1.0.10"

    conda "${projectDir}/environments/squirrel.yml"

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "squirrel/**", mode: 'copy', saveAs: {  fn -> fn.replace("squirrel/", "")}

    input:
        path fasta
        path refs
        path background
    output:
        path "squirrel/all_consensus.aln.fasta", emit: alignment
        path "squirrel/all_consensus.tree", emit: tree
        path "squirrel/**", emit: all
        path "squirrel.version", emit: version

  script:
    extra = ""
    if ( params.squirrel_assembly_refs )
        extra += " --assembly-refs ${refs}"
    if ( params.background_sequences )
        extra += " --background-file ${background}"
    if ( params.clade )
        extra += " --clade ${params.clade}"
    if ( params.run_phylo )
        extra += " --run-phylo --include-background"
    if ( params.run_apobec3_phylo )
        extra += " --run-apobec3-phylo --include-background"
    if ( params.outgroups )
        extra += " --outgroups ${params.outgroups}"

    """
    export XDG_CACHE_HOME=\$PWD/.cache
    squirrel --version 2>&1 | sed 's/: /,/' > squirrel.version
    squirrel ${fasta} --no-mask --seq-qc -o squirrel --outfile all_consensus.aln.fasta --tempdir squirrel_tmp -t ${task.cpus} ${extra}
    """
}
