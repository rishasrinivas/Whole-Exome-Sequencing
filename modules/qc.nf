// Quality Control modules

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.{zip,html}"), emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc -q -t ${task.cpus} ${reads}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path(files)

    output:
    path("multiqc_report.html"), emit: html
    path("multiqc_data"), emit: data

    script:
    """
    multiqc . \\
        --title "WES Pipeline Report" \\
        --filename multiqc_report.html \\
        --force
    """
}
