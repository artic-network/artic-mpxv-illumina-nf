process squirrelAlignmentAndQC {

    label "process_low"

    container "community.wave.seqera.io/library/squirrel:1.0.10--dc85d171b951f751"

    conda "bioconda::squirrel=1.0.10"

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
    if ( params.outgroups )
        extra += " --outgroups ${params.outgroups}"

    """
    export XDG_CACHE_HOME=\$PWD/.cache
    squirrel --version 2>&1 | sed 's/: /,/' > squirrel.version
    squirrel ${fasta} --no-mask --seq-qc -o squirrel --outfile all_consensus.aln.fasta --tempdir squirrel_tmp -t ${task.cpus} ${extra}
    """
}
