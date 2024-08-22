process performHostFilter {

    tag { sampleName }

    label 'process_medium'

    container 'community.wave.seqera.io/library/hostile:1.1.0--15df70fea624a735'

    conda 'bioconda::hostile=1.1.0'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}_hostfiltered_R*.fastq.gz", mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*.clean_1.fastq.gz"), path("*.clean_2.fastq.gz"), emit: fastqPairs

    script:
    """
    hostile clean --fastq1 ${forward} --fastq2 ${reverse} --out-dir . --threads ${task.cpus}
    """
}

process normalizeDepth {

    tag { sampleName }

    label  'process_medium'

    container 'community.wave.seqera.io/library/bbmap:39.06--5b971f29ed092959'

    conda 'bioconda::bbmap=39.06'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_norm_R{1,2}.fq.gz', mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("${sampleName}_norm_R1.fq.gz"), path("${sampleName}_norm_R2.fq.gz")

    script:
    """
    mkdir tmp
    memory=\$(echo '${task.memory}' | cut -d ' ' -f 1)

    bbnorm.sh \
      threads=${task.cpus} \
      target=${params.normalizationTargetDepth} \
      mindepth=${params.normalizationMinDepth} \
      in=${forward} \
      in2=${reverse} \
      out=${sampleName}_norm_R1.fq.gz \
      out2=${sampleName}_norm_R2.fq.gz \
      tmpdir=tmp \
    """
}

process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    label 'process_single'

    container 'community.wave.seqera.io/library/trim-galore:0.6.10--e1d78c153f940cdf'

    conda 'bioconda::trim-galore=0.6.10'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), optional: true

    script:
    """
    if [[ \$(gunzip -c ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      cp ${forward} ${sampleName}_hostfiltered_val_1.fq.gz
      cp ${reverse} ${sampleName}_hostfiltered_val_2.fq.gz
    else
      trim_galore --paired $forward $reverse
    fi
    """
}

process filterResidualAdapters {
    /**
    * Discard reads that contain residual adapter sequences that indicate trimming may have failed
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output untrim_filter_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    label 'process_single'

    container 'community.wave.seqera.io/library/pysam:0.22.1--6e6607269e942cf3'

    conda 'bioconda::pysam=0.22.1'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*{1,2}_posttrim_filter.fq.gz', mode: 'copy'

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*1_posttrim_filter.fq.gz"), path("*2_posttrim_filter.fq.gz")

    script:
    """
    filter_residual_adapters.py --input_R1 $forward --input_R2 $reverse
    """
}

process indexReference {
    /**
    * Indexes reference fasta file in the scheme repo using bwa.
    */

    tag { ref }

    label 'process_single'

    container 'community.wave.seqera.io/library/bwa:0.7.18--324359fbc6e00dba'

    conda 'bioconda::bwa=0.7.17'

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

    label 'process_medium'

    container 'community.wave.seqera.io/library/bwa_samtools:3938c84206f62975'

    conda 'bioconda::bwa=0.7.18', 'bioconda::samtools=1.12'

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

process trimPrimerSequences {

    tag { sampleName }

    label 'process_medium'

    container 'community.wave.seqera.io/library/ivar_samtools:e656be3eda151c7e'

    conda 'bioconda::ivar=1.4.2', 'bioconda::samtools=1.20'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped{.bam,.bam.bai}", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted{.bam,.bam.bai}", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bam_index), path(bedfile)

    output:
    tuple val(sampleName), path("${sampleName}.mapped.bam"), path("${sampleName}.mapped.bam.bai"), emit: mapped
    tuple val(sampleName), path("${sampleName}.mapped.primertrimmed.sorted.bam"), path("${sampleName}.mapped.primertrimmed.sorted.bam.bai" ), emit: ptrim

    script:
    """
    samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
    samtools index ${sampleName}.mapped.bam
    ivar trim -e -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.keepLen} -q ${params.qualThreshold} -p ivar.out
    samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
    samtools index ${sampleName}.mapped.primertrimmed.sorted.bam
    """
}

process callConsensusFreebayes {

    tag { sampleName }

    label 'process_single'

    container 'community.wave.seqera.io/library/bcftools_freebayes_pysam_tabix:c16cddc9a1a28b82'

    conda 'bioconda::freebayes=1.3.6', 'bioconda::bcftools=1.20', 'bioconda::pysam=0.22.1', 'bioconda::tabix=1.11'

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
    bcftools consensus -f ${ref} -I ${sampleName}.ambiguous.norm.vcf.gz > ${sampleName}.ambiguous.fa

    # Get viral contig name from reference
    CTG_NAME=\$(head -n1 ${ref} | sed 's/>//')

    # apply remaninng variants, including indels
    bcftools consensus -f ${sampleName}.ambiguous.fa -m ${sampleName}.mask.txt ${sampleName}.fixed.norm.vcf.gz | sed s/\$CTG_NAME/${sampleName}/ > ${sampleName}.consensus.fa
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

    conda 'bioconda::mafft=7.525'

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
