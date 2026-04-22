# @title RiboPause-NB: Local Z-score & Expression Filtering
# @description This module performs transcript-level expression filtering and calculates localized Z-scores 
# for ribosome stalling sites. By leveraging sliding-window statistics (mean and standard deviation) 
# across high-expression genes, it provides a normalized metric to identify significant pausing 
# events while suppressing noise and false positives in low-coverage regions.

# 1. Load Dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_packages <- c("Rsamtools", "GenomicAlignments", "IRanges", "data.table")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(IRanges)
  library(data.table)
})

#' @details Calculate local Z-scores and filter by gene expression levels
calculate_local_zscores <- function(psites, window_size = 100, expr_quantile = 0.5, fc_threshold = 10) {
  
  message("Step 1: Filtering by Gene Expression...")
  
  # Generate coverage track from P-sites
  cvg <- coverage(psites)
  
  # Calculate mean depth per transcript to identify high-expression genes
  gene_means <- mean(cvg)
  gene_expr <- data.table(gene = names(gene_means), mean_depth = as.numeric(gene_means))
  
  # Keep genes above the specified expression quantile (e.g., Top 50%)
  expr_threshold <- quantile(gene_expr[mean_depth > 0]$mean_depth, probs = expr_quantile)
  high_expr_genes <- gene_expr[mean_depth >= expr_threshold]$gene
  
  message(paste("   -> Keeping", length(high_expr_genes), "genes with mean depth >=", round(expr_threshold, 2)))
  
  # Filter coverage and counts to high-expression genes
  cvg_sub <- cvg[high_expr_genes]
  
  # Create a data table of observed counts at each occupied position
  dt_counts <- data.table(seqnames = as.character(seqnames(psites)), pos = start(psites))
  dt_counts <- dt_counts[seqnames %in% high_expr_genes, .(count = .N), by = .(seqnames, pos)]
  setorder(dt_counts, seqnames, pos)
  
  message("Step 2: Calculating Local Statistics using Sliding Windows...")
  
  # Define windows centered at each P-site
  target_gr <- GRanges(dt_counts$seqnames, IRanges(dt_counts$pos, width = 1))
  windows <- resize(target_gr, width = window_size, fix = "center")
  
  # Trim windows to stay within transcript boundaries
  # Note: seqlengths should be pre-assigned to psites from the reference FASTA
  windows <- trim(windows)
  
  # A. Calculate Local Mean (First Moment)
  v_x <- Views(cvg_sub, windows)
  sum_x <- viewSums(v_x)
  local_mean <- sum_x / window_size
  
  # B. Calculate Local Standard Deviation (Second Moment)
  # Var(X) = E[X^2] - (E[X])^2
  cvg_sq <- cvg_sub^2 
  v_x2 <- Views(cvg_sq, windows)
  sum_x2 <- viewSums(v_x2)
  mean_x2 <- sum_x2 / window_size
  
  local_var <- mean_x2 - (local_mean^2)
  local_sd <- sqrt(pmax(local_var, 0)) # pmax handles tiny negative values from float precision
  
  # --- Step 3: Compute Scores ---
  message("Step 3: Generating Final Z-scores and Fold Change...")
  
  dt_counts[, `:=`(
    window_mean = as.numeric(local_mean),
    window_sd   = as.numeric(local_sd)
  )]
  
  # 1. Fold Change (Pause Score) with regularization to avoid division by zero
  dt_counts[, pause_score := count / (window_mean + 1.0)]
  
  # 2. Local Z-score
  # Regularization (SD + 1.0) suppresses noise in low-variance regions
  dt_counts[, z_score := (count - window_mean) / (window_sd + 1.0)]
  
  # 3. Pause Ratio (Percentage of total window signal at this site)
  dt_counts[, pause_ratio_pct := (count / (as.numeric(sum_x) + 1.0)) * 100]
  
  # Filter based on significance thresholds
  final_res <- dt_counts[count > 0 & (pause_score >= fc_threshold | z_score >= 3.0)]
  
  message(paste("   -> Found", nrow(final_res), "significant pause sites."))
  
  return(final_res)
}
