# @title RiboPause-NB: P-site Processor Module
# @description This script serves as the foundational preprocessing module for the RiboPause-NB project. 
# It automates the extraction of P-site coordinates from Ribo-seq BAM files by integrating genomic 
# annotations and calculating length-specific offsets. The resulting dataset provides the essential 
# spatial features required for downstream Negative Binomial statistical modeling and machine 
# learning feature extraction.

# 1. Load Dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_packages <- c("riboWaltz", "GenomicFeatures", "data.table")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

suppressPackageStartupMessages({
  library(riboWaltz)
  library(GenomicFeatures)
  library(data.table)
})

#' @details Core P-site Processing Workflow
process_ribo_psites <- function(bam_dir, gtf_path, 
                                organism = "Mus musculus", 
                                dataSource = "Ensembl",
                                flanking = 6) {
  
  # Step 1: Generate Annotation Table 
  # Builds a transcriptome database from GTF for riboWaltz compatibility 
  message("Step 1: Creating annotation from GTF...")
  annotation_dt <- create_annotation(
    gtfpath = gtf_path,
    dataSource = dataSource,
    organism = organism
  )
  
  # Step 2: Convert BAM to List 
  # Loads BAM files and associated indices (.bai) from the specified folder 
  message("Step 2: Loading BAM files from directory...")
  reads_list <- bamtolist(
    bamfolder = bam_dir,
    annotation = annotation_dt,
    name_samples = NULL
  )
  
  # Step 3: Duplicate Filtering 
  # Removes reads sharing identical 5' and 3' coordinates to minimize PCR bias
  message("Step 3: Filtering duplicates...")
  filtered_list <- duplicates_filter(
    data = reads_list,
    extremity = "both"
  )
  
  # Step 4: P-site Offset Calculation
  # Identifies the ribosomal P-site position relative to the 5' or 3' end 
  # Filters out reads too close to the start codon based on flanking parameters
  message("Step 4: Calculating P-site offsets...")
  psite_offset <- psite(
    data = filtered_list,
    flanking = flanking,
    extremity = "auto"
  )
  
  # Step 5: Generate P-site Information
  # Applies calculated offsets to generate final P-site coordinates for the NB model
  message("Step 5: Generating P-site info...")
  reads_psite_list <- psite_info(
    data = filtered_list,
    offset = psite_offset
  )
  
  return(list(
    psite_data = reads_psite_list,
    annotation = annotation_dt,
    offsets = psite_offset
  ))
}
