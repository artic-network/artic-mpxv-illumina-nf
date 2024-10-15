process fetchHostileReference {

    label 'process_low'

    container 'community.wave.seqera.io/library/hostile:1.1.0--15df70fea624a735'

    conda 'bioconda::hostile=1.1.0'

    output:
    path "hostile.ok"

    script:
    """
    export HOSTILE_CACHE_DIR='${params.store_dir}/hostile'
    hostile fetch --aligner bowtie2

    touch hostile.ok
    """
}

process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    label 'process_low'

    container 'community.wave.seqera.io/library/trim-galore:0.6.10--e1d78c153f940cdf'

    conda 'bioconda::trim-galore=0.6.10'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*_val_1.fq.gz"), path("*_val_2.fq.gz")

    script:
    """
    trim_galore --cores ${task.cpus} --paired ${forward} ${reverse}
    """
}

process performHostFilter {

    tag { sampleName }

    label 'process_low'

    container 'community.wave.seqera.io/library/hostile:1.1.0--15df70fea624a735'

    conda 'bioconda::hostile=1.1.0'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.clean_*.fastq.gz", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)
    path hostile_ok

    output:
    tuple val(sampleName), path("*.clean_1.fastq.gz"), path("*.clean_2.fastq.gz")

    script:
    """
    export HOSTILE_CACHE_DIR='${params.store_dir}/hostile/'
    hostile clean --fastq1 ${forward} --fastq2 ${reverse} --out-dir . --threads ${task.cpus}
    """
}

process align_trim {

    tag { sampleName }

    label 'process_single'

    container 'community.wave.seqera.io/library/pysam_samtools_numpy_pandas:0ef969f9a905399f'

    conda 'bioconda::pysam=0.22.1 bioconda::samtools=1.12 conda-forge::numpy=2.1.1 conda-forge::pandas=2.2.2'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.rg.sorted.bam*", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.amplicon_depths.tsv", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bam_index)
    path bed

    output:
    tuple val(sampleName), path("${sampleName}.primertrimmed.rg.sorted.bam"), path("${sampleName}.primertrimmed.rg.sorted.bam.bai"), emit: ptrimmed_bam
    path "${sampleName}.amplicon_depths.tsv", emit: amplicon_depth_tsv

    script:
    if (!params.skip_normalize_depth) {
        normalise_string = "--normalise ${params.normalizationTargetDepth}"
    } else {
        normalise_string = ""
    }

    if (params.discard_incorrect_primer_pairs) {
        pp_string = "--discard-incorrect-primer-pairs"
    } else {
        pp_string = ""
    }

    """
    align_trim.py ${normalise_string} ${pp_string} ${bed} --paired --no-read-groups --primer-match-threshold ${params.primer_match_threshold} --min-mapq ${params.min_mapq} --trim-primers --report ${sampleName}.alignreport.csv --amp-depth-report ${sampleName}.amplicon_depths.tsv < ${bam} 2> ${sampleName}.alignreport.er | samtools sort -T ${sampleName} - -o ${sampleName}.primertrimmed.rg.sorted.bam && samtools index ${sampleName}.primertrimmed.rg.sorted.bam
    """
}

process indexReference {
    /**
    * Indexes reference fasta file in the scheme repo using bwa.
    */

    tag { ref }

    label 'process_single'

    container 'community.wave.seqera.io/library/bwa:0.7.18--324359fbc6e00dba'

    conda 'bioconda::bwa=0.7.18'

    input:
    path(ref)

    output:
    path "${ref}.*"

    script:
    """
    bwa index ${ref}
    """
}

process readMapping {
    /**
    * Maps trimmed paired fastq using BWA (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */

    tag { sampleName }

    label 'process_low'

    container 'community.wave.seqera.io/library/bwa_samtools:3938c84206f62975'

    conda 'bioconda::bwa=0.7.18 bioconda::samtools=1.12'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted{.bam,.bam.bai}", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse), path(ref)
    path bwaIndex

    output:
    tuple val(sampleName), path("${sampleName}.sorted.bam"), path("${sampleName}.sorted.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${forward} ${reverse} | \
    samtools sort -o ${sampleName}.sorted.bam
    samtools index ${sampleName}.sorted.bam
    """
}

process callConsensusFreebayes {

    tag { sampleName }

    label 'process_single'

    container 'community.wave.seqera.io/library/bcftools_freebayes_pysam_tabix:c16cddc9a1a28b82'

    conda 'bioconda::freebayes=1.3.6 bioconda::bcftools=1.20 bioconda::pysam=0.22.1 bioconda::tabix=1.11'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.consensus.fa", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.norm.vcf", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bam_index), path(ref)

    output:
    tuple val(sampleName), path("${sampleName}.consensus.fa"), emit: consensus
    tuple val(sampleName), path("${sampleName}.variants.norm.vcf"), emit: variants

    script:
    """
    # the sed is to fix the header until a release is made with this fix
    #   --gvcf --gvcf-dont-use-chunk true ${bam} | sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > ${sampleName}.gvcf
    # https://github.com/freebayes/freebayes/pull/549
    freebayes -p 1 \
              -f ${ref} \
              -F 0.2 \
              -C 1 \
              --pooled-continuous \
              --min-coverage ${params.varMinDepth} \
              --gvcf --gvcf-dont-use-chunk true ${bam} > ${sampleName}.gvcf
 
    # make depth mask, split variants into ambiguous/consensus
    # NB: this has to happen before bcftools norm or else the depth mask misses any bases exposed during normalization
    process_gvcf.py -d ${params.varMinDepth} \
                    -l ${params.varMinFreqThreshold} \
                    -u ${params.varFreqThreshold} \
                    -m ${sampleName}.mask.txt \
                    -v ${sampleName}.variants.vcf \
                    -c ${sampleName}.consensus.vcf ${sampleName}.gvcf

    # normalize variant records into canonical VCF representation
    for v in "variants" "consensus"; do
        bcftools norm -f ${ref} ${sampleName}.\$v.vcf > ${sampleName}.\$v.norm.vcf
    done

    # split the consensus sites file into a set that should be IUPAC codes and all other bases, using the ConsensusTag in the VCF
    for vt in "ambiguous" "fixed"; do
        cat ${sampleName}.consensus.norm.vcf | awk -v vartag=ConsensusTag=\$vt '\$0 ~ /^#/ || \$0 ~ vartag' > ${sampleName}.\$vt.norm.vcf
        bgzip -f ${sampleName}.\$vt.norm.vcf
        tabix -f -p vcf ${sampleName}.\$vt.norm.vcf.gz
    done

    # apply ambiguous variants first using IUPAC codes. this variant set cannot contain indels or the subsequent step will break
    bcftools consensus -s - -f ${ref} -I ${sampleName}.ambiguous.norm.vcf.gz > ${sampleName}.ambiguous.fa

    # Get viral contig name from reference
    CTG_NAME=\$(head -n1 ${ref} | sed 's/>//')

    # apply remaninng variants, including indels
    bcftools consensus -s - -f ${sampleName}.ambiguous.fa -m ${sampleName}.mask.txt ${sampleName}.fixed.norm.vcf.gz | sed s/\$CTG_NAME/${sampleName}/ > ${sampleName}.consensus.fa
    """
}

process alignConsensusToReference {
    /**
    * Aligns consensus sequence against reference using mafft. Uses the --keeplength
    * flag to guarantee that all alignments remain the same length as the reference.
    */

    tag { sampleName }

    label 'process_single'

    container 'community.wave.seqera.io/library/mafft:7.526--dbad4ff150905890'

    conda 'bioconda::mafft=7.526'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.consensus.aln.fa", mode: 'copy'

    input:
    tuple val(sampleName), path(consensus), path(reference)

    output:
    tuple val(sampleName), path("${sampleName}.consensus.aln.fa")

    script:
    // Convert multi-line fasta to single line
    awk_string = '/^>/ {printf("\\n%s\\n", $0); next; } { printf("%s", $0); }  END { printf("\\n"); }'
    """
    mafft \
      --preservecase \
      --keeplength \
      --add \
      ${consensus} \
      ${reference} \
      > ${sampleName}.with_ref.multi_line.alignment.fa
    awk '${awk_string}' ${sampleName}.with_ref.multi_line.alignment.fa > ${sampleName}.with_ref.single_line.alignment.fa
    tail -n 2 ${sampleName}.with_ref.single_line.alignment.fa > ${sampleName}.consensus.aln.fa
    """
}

process annotateVariantsVCF {
    /**
    */

    tag { sampleName }

    label 'process_single'

    container 'biocontainers/bcftools:1.20--h8b25389_1'

    conda 'bioconda::bcftools=1.20'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.norm.consequence.{vcf,tsv}", mode: 'copy'

    input:
        tuple val(sampleName), path(vcf), path(ref), path(gff)

    output:
        tuple val(sampleName), path("${sampleName}.variants.norm.consequence.vcf"), emit: vcf
        tuple val(sampleName), path("${sampleName}.variants.norm.consequence.tsv"), emit: tsv

    script:
        """
        bcftools csq -f ${ref} -g ${gff} ${vcf} -Ov -o ${sampleName}.variants.norm.consequence.vcf
        bcftools csq -f ${ref} -g ${gff} ${vcf} -Ot -o ${sampleName}.variants.norm.consequence.tsv
        """
}
