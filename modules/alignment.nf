// Alignment and BAM processing modules

import java.nio.file.Paths

process BWA_MEM2_INDEX {
    publishDir "${params.outdir}/reference", mode: 'copy', saveAs: { filename -> 
        params.save_intermediate ? filename : null 
    }

    input:
    path(fasta)

    output:
    tuple path(fasta), path("${fasta}.*"), emit: index

    script:
    """
    bwa-mem2 index ${fasta}
    """
}

process BWA_MEM2_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.bam*", saveAs: { filename -> 
        params.save_intermediate ? filename : null 
    }

    input:
    tuple val(sample_id), path(reads)
    tuple path(ref_fasta), path(ref_index)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), emit: bam

    script:
    def read_group = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:WES\\tPL:ILLUMINA\\tPU:unit1"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${read_group}" \\
        ${ref_fasta} \\
        ${reads[0]} \\
        ${reads[1]} | \\
    samtools view -@ ${task.cpus} -bS - > ${sample_id}.aligned.bam
    """
}

process PICARD_SORT {
    tag "$sample_id"
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.sorted.ba*", saveAs: { filename -> 
        params.save_intermediate ? filename : null 
    }

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bai"), emit: bam

    script:
    """
    picard SortSam \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}.sorted.bam \\
        SORT_ORDER=coordinate \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        MAX_RECORDS_IN_RAM=2000000
    """
}

process PICARD_MARKDUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bai"), emit: bam
    path("${sample_id}.dedup.metrics.txt"), emit: metrics

    script:
    """
    picard MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}.dedup.bam \\
        METRICS_FILE=${sample_id}.dedup.metrics.txt \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        ASSUME_SORT_ORDER=coordinate \\
        MAX_RECORDS_IN_RAM=2000000
    """
}
