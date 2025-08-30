// Report generation modules

process GENERATE_REPORT {
    tag "$sample_id"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(sample_id), path(annotated_vcf)
    path(fastqc_files)
    path(dedup_metrics)

    output:
    path("${sample_id}_clinical_report.html"), emit: html_report
    path("${sample_id}_clinical_report.txt"), emit: text_report
    path("${sample_id}_summary.json"), emit: json_summary
    path("reports/"), emit: reports_dir

    script:
    """
    mkdir -p reports
    
    python3 ${projectDir}/generate_report.py \\
        ${annotated_vcf} \\
        --output reports \\
        --sample-name ${sample_id} \\
        --format all
    
    # Copy reports to main directory for publishing
    cp reports/* ./ || true
    """
}
