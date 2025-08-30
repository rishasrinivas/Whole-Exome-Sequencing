#!/usr/bin/env nextflow

/*
========================================================================================
    WES (Whole Exome Sequencing) Analysis Pipeline
========================================================================================
    Complete pipeline for WES data analysis from FASTQ to clinical report
    
    Author: WES Pipeline Team
    Version: 1.0.0
========================================================================================
*/

nextflow.enable.dsl = 2

// Print help message
def helpMessage() {
    log.info"""
    ========================================
     WES Analysis Pipeline v${workflow.manifest.version}
    ========================================
    
    Usage:
    nextflow run main.nf --input_dir <path> --ref_genome_fasta <path> [options]
    
    Required Arguments:
      --input_dir                   Path to directory containing FASTQ files
      --ref_genome_fasta           Path to reference genome FASTA file
      
    Optional Arguments:
      --outdir                     Output directory [default: ./results]
      --target_bed                 BED file with target regions for exome
      --dbsnp                      dbSNP VCF file for variant calling
      --known_indels               Known indels VCF file for GATK
      --vep_cache_dir             VEP cache directory
      --sample_sheet              Sample sheet CSV file
      --help                       Show this help message
    
    Example:
    nextflow run main.nf \\
        --input_dir ./data/raw \\
        --ref_genome_fasta ./reference/GRCh38.fa \\
        --target_bed ./reference/exome_targets.bed \\
        --outdir ./results
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.ref_genome_fasta) {
    error "Error: --ref_genome_fasta is required"
}

if (!params.input_dir) {
    error "Error: --input_dir is required"
}

// Import modules
include { validateInputs } from './modules/utils'
include { FASTQC } from './modules/qc'
include { MULTIQC } from './modules/qc'
include { BWA_MEM2_INDEX } from './modules/alignment'
include { BWA_MEM2_ALIGN } from './modules/alignment'
include { PICARD_SORT } from './modules/alignment'
include { PICARD_MARKDUPLICATES } from './modules/alignment'
include { GATK_CREATESEQUENCEDICTIONARY } from './modules/preprocessing'
include { SAMTOOLS_FAIDX } from './modules/preprocessing'
include { GATK_HAPLOTYPECALLER } from './modules/variant_calling'
include { VEP_ANNOTATE } from './modules/annotation'
include { GENERATE_REPORT } from './modules/reporting'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    // Validate inputs
    validateInputs()
    
    // Create input channel from FASTQ files
    fastq_ch = Channel
        .fromFilePairs("${params.input_dir}/*_{R1,R2,1,2}*.{fastq,fq}.gz", checkIfExists: true)
        .ifEmpty { error "Cannot find any FASTQ files in ${params.input_dir}" }
    
    // Reference genome files
    ref_fasta = file(params.ref_genome_fasta, checkIfExists: true)
    
    // Prepare reference genome
    SAMTOOLS_FAIDX(ref_fasta)
    GATK_CREATESEQUENCEDICTIONARY(ref_fasta)
    BWA_MEM2_INDEX(ref_fasta)
    
    // Quality control
    FASTQC(fastq_ch)
    
    // Alignment
    BWA_MEM2_ALIGN(
        fastq_ch,
        BWA_MEM2_INDEX.out.index
    )
    
    // Sort BAM files
    PICARD_SORT(BWA_MEM2_ALIGN.out.bam)
    
    // Mark duplicates
    PICARD_MARKDUPLICATES(PICARD_SORT.out.bam)
    
    // Variant calling
    GATK_HAPLOTYPECALLER(
        PICARD_MARKDUPLICATES.out.bam,
        ref_fasta,
        SAMTOOLS_FAIDX.out.fai,
        GATK_CREATESEQUENCEDICTIONARY.out.dict
    )
    
    // Variant annotation
    VEP_ANNOTATE(GATK_HAPLOTYPECALLER.out.vcf)
    
    // Generate final report
    GENERATE_REPORT(
        VEP_ANNOTATE.out.annotated_vcf,
        FASTQC.out.zip.collect(),
        PICARD_MARKDUPLICATES.out.metrics.collect()
    )
    
    // MultiQC report
    if (!params.skip_multiqc) {
        multiqc_input = Channel.empty()
            .mix(FASTQC.out.zip.collect())
            .mix(PICARD_MARKDUPLICATES.out.metrics.collect())
        
        MULTIQC(multiqc_input.collect())
    }
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info"""
    ========================================
     Pipeline Completed Successfully!
    ========================================
    
    Output directory: ${params.outdir}
    
    Key outputs:
    - QC Reports: ${params.outdir}/qc/
    - Aligned BAMs: ${params.outdir}/alignment/
    - Variant Calls: ${params.outdir}/variants/
    - Annotations: ${params.outdir}/annotation/
    - Final Report: ${params.outdir}/reports/
    
    Pipeline execution summary:
    - Completed at: ${workflow.complete}
    - Duration: ${workflow.duration}
    - Success: ${workflow.success}
    - Work directory: ${workflow.workDir}
    """.stripIndent()
}

workflow.onError {
    log.error"""
    ========================================
     Pipeline Execution Failed
    ========================================
    
    Error message: ${workflow.errorMessage}
    Failed process: ${workflow.errorReport}
    Exit status: ${workflow.exitStatus}
    Work directory: ${workflow.workDir}
    """.stripIndent()
}
