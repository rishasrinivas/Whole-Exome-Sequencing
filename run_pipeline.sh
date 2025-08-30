#!/bin/bash

#======================================
# WES Pipeline Launcher Script
#======================================

set -e

# Default parameters
INPUT_DIR=""
REF_GENOME=""
OUTPUT_DIR="./results"
PROFILE="conda"
RESUME=""
TARGET_BED=""
DBSNP=""
VEP_CACHE=""
MAX_CPUS=8
MAX_MEMORY="32.GB"

# Help function
show_help() {
    cat << EOF
WES Pipeline Launcher

Usage: $0 [OPTIONS]

Required Options:
  -i, --input-dir PATH        Input directory containing FASTQ files
  -r, --reference PATH        Reference genome FASTA file

Optional Options:
  -o, --output-dir PATH       Output directory [default: ./results]
  -p, --profile NAME          Execution profile [default: conda]
  -t, --target-bed PATH       Target regions BED file
  -d, --dbsnp PATH           dbSNP VCF file
  -c, --vep-cache PATH       VEP cache directory
  --max-cpus INT             Maximum CPU cores [default: 8]
  --max-memory SIZE          Maximum memory [default: 32.GB]
  --resume                   Resume previous run
  -h, --help                 Show this help message

Profiles:
  conda                      Use conda environment
  docker                     Use Docker containers
  singularity               Use Singularity containers
  cluster                    SLURM cluster execution
  test                       Run with test data

Examples:
  # Basic usage
  $0 -i ./data/raw -r ./reference/GRCh38.fa
  
  # Full analysis with all options
  $0 -i ./data/raw \\
     -r ./reference/GRCh38.fa \\
     -t ./reference/exome_targets.bed \\
     -d ./reference/dbsnp.vcf.gz \\
     -c ./vep_cache \\
     -o ./my_results \\
     --max-cpus 16 \\
     --max-memory 64.GB
  
  # Resume interrupted run
  $0 -i ./data/raw -r ./reference/GRCh38.fa --resume
  
  # Test run
  $0 -p test

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -r|--reference)
            REF_GENOME="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -t|--target-bed)
            TARGET_BED="$2"
            shift 2
            ;;
        -d|--dbsnp)
            DBSNP="$2"
            shift 2
            ;;
        -c|--vep-cache)
            VEP_CACHE="$2"
            shift 2
            ;;
        --max-cpus)
            MAX_CPUS="$2"
            shift 2
            ;;
        --max-memory)
            MAX_MEMORY="$2"
            shift 2
            ;;
        --resume)
            RESUME="-resume"
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Check for test profile
if [[ "$PROFILE" == "test" ]]; then
    echo "Running test profile..."
    nextflow run main.nf -profile test $RESUME
    exit 0
fi

# Validate required parameters
if [[ -z "$INPUT_DIR" ]] && [[ "$PROFILE" != "test" ]]; then
    echo "Error: Input directory is required (-i/--input-dir)"
    show_help
    exit 1
fi

if [[ -z "$REF_GENOME" ]] && [[ "$PROFILE" != "test" ]]; then
    echo "Error: Reference genome is required (-r/--reference)"
    show_help
    exit 1
fi

# Check if files exist
if [[ -n "$INPUT_DIR" ]] && [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [[ -n "$REF_GENOME" ]] && [[ ! -f "$REF_GENOME" ]]; then
    echo "Error: Reference genome not found: $REF_GENOME"
    exit 1
fi

if [[ -n "$TARGET_BED" ]] && [[ ! -f "$TARGET_BED" ]]; then
    echo "Error: Target BED file not found: $TARGET_BED"
    exit 1
fi

if [[ -n "$DBSNP" ]] && [[ ! -f "$DBSNP" ]]; then
    echo "Error: dbSNP file not found: $DBSNP"
    exit 1
fi

if [[ -n "$VEP_CACHE" ]] && [[ ! -d "$VEP_CACHE" ]]; then
    echo "Error: VEP cache directory not found: $VEP_CACHE"
    exit 1
fi

# Build command
CMD="nextflow run main.nf"
CMD="$CMD --input_dir $INPUT_DIR"
CMD="$CMD --ref_genome_fasta $REF_GENOME"
CMD="$CMD --outdir $OUTPUT_DIR"
CMD="$CMD --max_cpus $MAX_CPUS"
CMD="$CMD --max_memory $MAX_MEMORY"
CMD="$CMD -profile $PROFILE"

# Add optional parameters
if [[ -n "$TARGET_BED" ]]; then
    CMD="$CMD --target_bed $TARGET_BED"
fi

if [[ -n "$DBSNP" ]]; then
    CMD="$CMD --dbsnp $DBSNP"
fi

if [[ -n "$VEP_CACHE" ]]; then
    CMD="$CMD --vep_cache_dir $VEP_CACHE"
fi

# Add resume if specified
if [[ -n "$RESUME" ]]; then
    CMD="$CMD $RESUME"
fi

# Print configuration
echo "======================================="
echo "WES Pipeline Configuration"
echo "======================================="
echo "Input Directory: $INPUT_DIR"
echo "Reference Genome: $REF_GENOME"
echo "Output Directory: $OUTPUT_DIR"
echo "Profile: $PROFILE"
echo "Max CPUs: $MAX_CPUS"
echo "Max Memory: $MAX_MEMORY"
if [[ -n "$TARGET_BED" ]]; then
    echo "Target BED: $TARGET_BED"
fi
if [[ -n "$DBSNP" ]]; then
    echo "dbSNP: $DBSNP"
fi
if [[ -n "$VEP_CACHE" ]]; then
    echo "VEP Cache: $VEP_CACHE"
fi
if [[ -n "$RESUME" ]]; then
    echo "Resume: Yes"
fi
echo "======================================="
echo

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the pipeline
echo "Starting pipeline execution..."
echo "Command: $CMD"
echo

eval $CMD

echo
echo "Pipeline execution completed!"
echo "Results available in: $OUTPUT_DIR"
