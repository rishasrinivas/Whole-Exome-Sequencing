// Utility functions for WES pipeline

def validateInputs() {
    // Check if required parameters are provided
    if (!params.ref_genome_fasta) {
        error "Reference genome FASTA file is required (--ref_genome_fasta)"
    }
    
    if (!file(params.ref_genome_fasta).exists()) {
        error "Reference genome FASTA file not found: ${params.ref_genome_fasta}"
    }
    
    if (!params.input_dir) {
        error "Input directory is required (--input_dir)"
    }
    
    if (!file(params.input_dir).exists()) {
        error "Input directory not found: ${params.input_dir}"
    }
    
    // Check optional files
    if (params.target_bed && !file(params.target_bed).exists()) {
        error "Target BED file not found: ${params.target_bed}"
    }
    
    if (params.dbsnp && !file(params.dbsnp).exists()) {
        error "dbSNP file not found: ${params.dbsnp}"
    }
    
    if (params.known_indels && !file(params.known_indels).exists()) {
        error "Known indels file not found: ${params.known_indels}"
    }
    
    if (params.vep_cache_dir && !file(params.vep_cache_dir).exists()) {
        error "VEP cache directory not found: ${params.vep_cache_dir}"
    }
    
    log.info "Input validation passed"
}

def getSampleName(read_pair) {
    // Extract sample name from file name
    def (sample_id, path) = read_pair
    return sample_id
}

def getReadGroup(sample_name, read_pair) {
    // Generate read group string for BWA
    def rg_id = sample_name
    def rg_sm = sample_name
    def rg_lb = "WES"
    def rg_pl = "ILLUMINA"
    def rg_pu = "unit1"
    
    return "@RG\\tID:${rg_id}\\tSM:${rg_sm}\\tLB:${rg_lb}\\tPL:${rg_pl}\\tPU:${rg_pu}"
}
