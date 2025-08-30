# WES (Whole Exome Sequencing) Analysis Pipeline

A complete, production-ready Nextflow pipeline for whole exome sequencing data analysis, from raw FASTQ files to clinically-relevant reports.

## Overview

This pipeline performs comprehensive WES analysis including:

1. **Quality Control** - FastQC and MultiQC reports
2. **Read Alignment** - BWA-MEM2 alignment to reference genome
3. **Variant Calling** - GATK HaplotypeCaller for SNVs and indels
4. **Variant Annotation** - Ensembl VEP with ClinVar, dbSNP, and gnomAD
5. **Clinical Reporting** - Automated generation of clinical reports highlighting pathogenic variants

## Features

- ✅ **Production-ready**: Robust error handling and resource management
- ✅ **Scalable**: Configurable for local, cluster, or cloud execution
- ✅ **Comprehensive**: Full pipeline from FASTQ to clinical report
- ✅ **Clinical focus**: Automated filtering and reporting of pathogenic variants
- ✅ **Quality controlled**: Integrated QC at every step
- ✅ **Reproducible**: Containerized environment with version control

## Quick Start

### Prerequisites

- Nextflow ≥23.10.0
- Conda or Docker/Singularity
- Minimum 16GB RAM, 8 CPU cores recommended
- ~100GB storage for reference data and outputs

### Installation

1. **Clone the repository**
```bash
git clone <repository-url>
cd wes-pipeline
```

2. **Set up the environment**
```bash
conda env create -f conda.yml
conda activate wes_pipeline
```

3. **Download reference data**
```bash
# Download human reference genome (GRCh38)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna reference/GRCh38.fa

# Download dbSNP (optional but recommended)
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz
mv GCF_000001405.39.gz reference/dbsnp.vcf.gz

# Download exome target regions
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/exome_calling_regions.v1.interval_list
# Convert to BED format if needed
```

4. **Prepare VEP cache (recommended)**
```bash
# Download VEP cache for better performance
vep_install -a cf -s homo_sapiens -y GRCh38 -c vep_cache/
```

### Basic Usage

**Minimal command:**
```bash
nextflow run main.nf \
    --input_dir ./data/raw \
    --ref_genome_fasta ./reference/GRCh38.fa \
    --outdir ./results
```

**Full command with all options:**
```bash
nextflow run main.nf \
    --input_dir ./data/raw \
    --ref_genome_fasta ./reference/GRCh38.fa \
    --target_bed ./reference/exome_targets.bed \
    --dbsnp ./reference/dbsnp.vcf.gz \
    --vep_cache_dir ./vep_cache \
    --outdir ./results \
    --max_cpus 16 \
    --max_memory '32.GB' \
    -profile conda
```

## Input Requirements

### Required Inputs

1. **FASTQ files**: Paired-end sequencing files
   - Format: `*_{R1,R2,1,2}*.{fastq,fq}.gz`
   - Location: Specified by `--input_dir`
   - Example: `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`

2. **Reference genome**: Human reference genome in FASTA format
   - Recommended: GRCh38/hg38
   - Must be uncompressed
   - Specify with `--ref_genome_fasta`

### Optional Inputs

1. **Target regions BED file** (`--target_bed`)
   - Exome capture regions
   - Improves variant calling accuracy and speed

2. **dbSNP VCF file** (`--dbsnp`)
   - Known variant annotations
   - Improves variant calling quality

3. **VEP cache directory** (`--vep_cache_dir`)
   - Pre-downloaded VEP annotations
   - Faster than online queries

## Directory Structure

```
wes-pipeline/
├── main.nf                 # Main pipeline script
├── nextflow.config         # Pipeline configuration
├── conda.yml              # Software dependencies
├── generate_report.py      # Clinical report generator
├── modules/                # Pipeline modules
│   ├── utils.nf
│   ├── qc.nf
│   ├── preprocessing.nf
│   ├── alignment.nf
│   ├── variant_calling.nf
│   ├── annotation.nf
│   └── reporting.nf
├── data/                   # Input data directory
│   └── raw/                # Raw FASTQ files
├── reference/              # Reference files
│   ├── GRCh38.fa
│   ├── dbsnp.vcf.gz
│   └── exome_targets.bed
└── results/                # Output directory
    ├── qc/                 # Quality control reports
    ├── alignment/          # Aligned BAM files
    ├── variants/           # Variant calls (VCF)
    ├── annotation/         # Annotated variants
    ├── reports/            # Clinical reports
    └── pipeline_info/      # Pipeline execution info
```

## Output Files

### Quality Control
- `results/qc/fastqc/`: FastQC reports for each sample
- `results/qc/multiqc/`: Aggregated quality report

### Alignment
- `results/alignment/`: Sorted, deduplicated BAM files
- `*.dedup.metrics.txt`: Duplication statistics

### Variant Calling
- `results/variants/`: Raw variant calls in VCF format
- `*.vcf.gz` and `*.vcf.gz.tbi`: Compressed VCF with index

### Annotation
- `results/annotation/`: VEP-annotated variants
- `*.annotated.vcf.gz`: Annotated VCF files
- `*.vep_summary.html`: VEP annotation summary

### Clinical Reports
- `results/reports/`: Final clinical reports
- `*_clinical_report.html`: Interactive HTML report
- `*_clinical_report.txt`: Text-based summary
- `*_summary.json`: Machine-readable summary

## Configuration Options

### Resource Configuration

```bash
# Adjust resources based on your system
--max_cpus 32              # Maximum CPU cores
--max_memory '128.GB'      # Maximum memory
--max_time '240.h'         # Maximum runtime
```

### Quality Filters

```bash
--min_dp 10               # Minimum read depth
--min_gq 20              # Minimum genotype quality
--max_missing 0.1        # Maximum missing data
```

### Pipeline Behavior

```bash
--save_intermediate      # Save intermediate files
--skip_multiqc          # Skip MultiQC report generation
```

## Execution Profiles

### Local Execution
```bash
nextflow run main.nf -profile standard
```

### With Conda
```bash
nextflow run main.nf -profile conda
```

### With Docker
```bash
nextflow run main.nf -profile docker
```

### Cluster Execution (SLURM)
```bash
nextflow run main.nf -profile cluster
```

### Test Run
```bash
nextflow run main.nf -profile test
```

## Clinical Report Interpretation

The pipeline generates comprehensive clinical reports focusing on:

### Key Sections

1. **Executive Summary**
   - Total variant counts
   - Pathogenic variant summary
   - Quality metrics overview

2. **Clinically Significant Variants**
   - Pathogenic and likely pathogenic variants
   - Gene information and protein changes
   - Variant quality metrics

3. **Methodology**
   - Analysis pipeline details
   - Quality thresholds
   - Reference databases used

### Variant Classifications

- **Pathogenic**: Disease-causing variants requiring clinical action
- **Likely Pathogenic**: Probably disease-causing, clinical correlation advised  
- **Uncertain Significance**: Unknown clinical relevance
- **Likely Benign/Benign**: Probably not disease-causing

### Report Formats

1. **HTML Report**: Interactive, web-viewable format
2. **Text Report**: Simple text format for integration
3. **JSON Summary**: Machine-readable data for downstream analysis

## Troubleshooting

### Common Issues

**1. Out of Memory Errors**
```bash
# Increase memory allocation
--max_memory '64.GB'
--gatk_java_opts '-Xmx8g'
```

**2. VEP Annotation Fails**
```bash
# Use offline cache for reliability
--vep_cache_dir ./vep_cache
```

**3. No Variants Found**
- Check target BED file format
- Verify reference genome version
- Review quality filters

**4. Pipeline Fails to Resume**
```bash
# Clean work directory
rm -rf work/
```

### Log Files

Pipeline execution information:
- `results/pipeline_info/execution_report.html`
- `results/pipeline_info/execution_timeline.html`
- `results/pipeline_info/execution_trace.txt`

### Support

For issues and questions:
1. Check the log files for error messages
2. Verify input file formats and paths
3. Ensure sufficient system resources
4. Review parameter configurations

## Performance Guidelines

### System Requirements

| Sample Count | CPU Cores | Memory | Storage | Runtime |
|-------------|-----------|--------|---------|---------|
| 1 sample    | 8 cores   | 32GB   | 50GB    | 4-8h    |
| 5 samples   | 16 cores  | 64GB   | 200GB   | 8-16h   |
| 10 samples  | 32 cores  | 128GB  | 500GB   | 12-24h  |

### Optimization Tips

1. **Use SSD storage** for work directory
2. **Pre-download VEP cache** for faster annotation  
3. **Use target BED file** to focus on exonic regions
4. **Adjust resource allocation** based on sample count
5. **Enable resume** for interrupted runs (`-resume`)

## Pipeline Validation

### Test Dataset

A small test dataset is included for validation:
```bash
nextflow run main.nf -profile test
```

### Expected Outputs

Successful completion should generate:
- QC reports with PASS status
- Aligned BAM files with >95% mapping rate
- Variant calls with appropriate quality scores
- Clinical report with variant annotations

## Advanced Usage

### Custom Annotation

Add custom VEP plugins:
```groovy
// In nextflow.config
params.vep_plugins = [
    'CADD': '/path/to/CADD',
    'dbNSFP': '/path/to/dbNSFP'
]
```

### Sample Sheet Input

Use CSV sample sheet for complex datasets:
```csv
sample_id,read1,read2,phenotype
Sample1,/path/to/R1.fq.gz,/path/to/R2.fq.gz,affected
Sample2,/path/to/R1.fq.gz,/path/to/R2.fq.gz,control
```

### Integration with Other Pipelines

The pipeline outputs standard formats (BAM, VCF) compatible with:
- GATK workflows
- Variant annotation pipelines  
- Clinical reporting systems
- Population genetics tools

## License

This pipeline is released under the MIT License. See LICENSE file for details.

## Citations

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319 (2017).
- **BWA-MEM2**: Vasimuddin, M. et al. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (2019).
- **GATK**: McKenna, A. et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20, 1297–1303 (2010).
- **VEP**: McLaren, W. et al. The Ensembl Variant Effect Predictor. Genome Biol. 17, 122 (2016).

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## Changelog

### Version 1.0.0
- Initial release
- Complete WES analysis pipeline
- Clinical report generation
- Multi-profile execution support
