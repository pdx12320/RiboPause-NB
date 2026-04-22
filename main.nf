// Define default parameters
params.bam_files = "data/*.bam"
params.gtf = "data/annotation.gtf"
params.outdir = "results"

process PROCESS_PSITES {
    publishDir "${params.outdir}/psites", mode: 'copy'
    
    input:
    path bam
    path gtf

    output:
    path "${bam.baseName}_psite.rds"

    script:
    """
    Rscript p_site_processor.R --bam $bam --annot $gtf --output ${bam.baseName}_psite.rds
    """
}

process RUN_NB_MODEL {
    publishDir "${params.outdir}/stalling", mode: 'copy'

    input:
    path psite_data

    output:
    path "${psite_data.baseName}_stalling.csv"

    script:
    """
    Rscript nb_model.R --input $psite_data --output ${psite_data.baseName}_stalling.csv
    """
}

workflow {
    // Create a channel from all BAM files found in the path
    bam_ch = Channel.fromPath(params.bam_files)
    
    // Workflow logic: Each BAM flows through the processes independently
    psite_ch = PROCESS_PSITES(bam_ch, params.gtf)
    RUN_NB_MODEL(psite_ch)
}
