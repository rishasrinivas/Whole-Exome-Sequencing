// Variant calling modules

process GATK_HAPLOTYPECALLER {
    tag "$sample_id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(ref_fasta)
    path(ref_fai)
    path(ref_dict)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf

    script:
    def intervals_arg = params.target_bed ? "-L ${params.target_bed}" : ""
    def dbsnp_arg = params.dbsnp ? "--dbsnp ${params.dbsnp}" : ""
    def java_opts = params.gatk_java_opts ?: "-Xmx4g"
    """
    gatk --java-options "${java_opts}" HaplotypeCaller \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        -O ${sample_id}.vcf.gz \\
        ${intervals_arg} \\
        ${dbsnp_arg} \\
        --sample-ploidy ${params.ploidy} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --create-output-variant-index
    """
}
