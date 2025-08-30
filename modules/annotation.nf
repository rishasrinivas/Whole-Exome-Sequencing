// Variant annotation modules

process VEP_ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf.gz"), emit: annotated_vcf
    path("${sample_id}.vep_summary.html"), emit: summary

    script:
    def cache_arg = params.vep_cache_dir ? "--cache --dir_cache ${params.vep_cache_dir}" : "--database"
    def species = params.vep_species ?: "homo_sapiens"
    def assembly = params.vep_assembly ?: "GRCh38"
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}.annotated.vcf.gz \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --species ${species} \\
        --assembly ${assembly} \\
        ${cache_arg} \\
        --fork ${task.cpus} \\
        --everything \\
        --check_existing \\
        --pick \\
        --plugin CADD,/path/to/CADD \\
        --custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \\
        --stats_file ${sample_id}.vep_summary.html \\
        --warning_file ${sample_id}.vep_warnings.txt \\
        --no_progress
    
    # Index the output VCF
    tabix -p vcf ${sample_id}.annotated.vcf.gz
    """
}
