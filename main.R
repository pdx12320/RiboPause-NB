# @title RiboPause-NB: Main Analysis Pipeline
# @description This is the primary entry point for the RiboPause-NB pipeline. It orchestrates the entire 
# workflow, from raw BAM preprocessing and P-site mapping to advanced statistical modeling (Negative Binomial) 
# and ML-ready feature extraction. The pipeline integrates local Z-score normalization with global 
# dispersion fitting to provide a comprehensive landscape of ribosome stalling events.

# 1. Load Essential Orchestration Packages
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# 2. Source Individual Modules
# Assuming all module scripts are in the same directory
source("p_site_processor.R")
source("z-score.R")
source("nb_model.R")
source("context_extraction.R")

# 3. Command Line Argument Parsing
option_list <- list(
  make_option(c("-b", "--bam"), type="character", help="Path to folder containing BAM and .bai files"),
  make_option(c("-r", "--ref"), type="character", help="Path to reference cDNA FASTA file"),
  make_option(c("-a", "--annot"), type="character", help="Path to GTF annotation file"),
  make_option(c("-l", "--lengths"), type="character", default="28,29,30", help="Comma-separated read lengths to analyze"),
  make_option(c("-o", "--offset"), type="character", default="12,12,12", help="Comma-separated P-site offsets corresponding to lengths"),
  make_option(c("-w", "--window"), type="integer", default=100, help="Sliding window size for local statistics"),
  make_option(c("-f", "--foldchange"), type="numeric", default=10.0, help="Fold-change threshold for stalling detection"),
  make_option(c("-q", "--quantile"), type="numeric", default=0.5, help="Expression quantile for gene filtering (0-1)"),
  make_option(c("-p", "--output"), type="character", default="stalling_results.csv", help="Output filename")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Basic Validation
if (is.null(opt$bam) || is.null(opt$ref) || is.null(opt$annot)) {
  stop("BAM directory, Reference FASTA, and GTF Annotation are required. Use --help for usage.")
}

# 4. Execution Workflow
main <- function() {
  
  message("=== Starting RiboPause-NB Pipeline ===")
  
  # Step 1: P-site Processing
  # Converts BAM to P-site coordinates using riboWaltz logic
  psite_results <- process_ribo_psites(
    bam_dir = opt$bam, 
    gtf_path = opt$annot, 
    organism = "Mus musculus" # Can be parameterized if needed
  )
  
  # Extract P-sites for the first sample (or combine if multi-sample)
  # For simplicity, this template processes the first sample in the list
  sample_name <- names(psite_results$psite_data)[1]
  psites_gr <- psite_results$psite_data[[sample_name]]
  
  # Step 2: Expression Filtering & Local Z-scores
  # Identifies high-expression genes and computes window-based metrics
  final_res <- calculate_local_zscores(
    psites = psites_gr, 
    window_size = opt$window, 
    expr_quantile = opt$quantile, 
    fc_threshold = opt$foldchange
  )
  
  # Step 3: Negative Binomial Modeling
  # Fits dispersion and calculates statistical significance (P-values/FDR)
  final_res <- fit_nb_stalling_model(
    final_res = final_res, 
    z_window = opt$window
  )
  
  # Step 4: Context & ML Feature Extraction
  # Fetches sequences and generates padded score tracks for downstream ML
  final_res <- extract_context_features(
    final_res = final_res, 
    fa_file = opt$ref, 
    psite_counts = psites_gr, # Use raw coverage for contextual score tracks
    up_len = 50, 
    down_len = 50
  )
  
  # Step 5: Final In-frame Annotation & Save
  # Map to codons and save the final ML-ready dataset
  message("Step 5: Finalizing and Exporting Results...")
  
  # Optional: Final filtering for Significant & In-frame sites
  # (Custom logic can be added here based on padj or z_score)
  
  write.csv(final_res, file = opt$output, row.names = FALSE, quote = TRUE)
  
  message(paste("=== Pipeline Complete! Results saved to:", opt$output, "==="))
}

# Run the pipeline
main()
