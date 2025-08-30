// Reference genome preprocessing modules

process SAMTOOLS_FAIDX {
    publishDir "${params.outdir}/reference", mode: 'copy', saveAs: { filename -> 
        params.save_intermediate ? filename : null 
    }

    input:
    path(fasta)

    output:
    path("${fasta}.fai"), emit: fai

    script:
    """
    samtools faidx ${fasta}
    """
}

process GATK_CREATESEQUENCEDICTIONARY {
    publishDir "${params.outdir}/reference", mode: 'copy', saveAs: { filename -> 
        params.save_intermediate ? filename : null 
    }

    input:
    path(fasta)

    output:
    path("${fasta.baseName}.dict"), emit: dict

    script:
    """
    gatk CreateSequenceDictionary \\
        -R ${fasta} \\
        -O ${fasta.baseName}.dict
    """
}
