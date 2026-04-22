# @title RiboPause-NB: Negative Binomial Modeling & Statistical Inference
# @description This module implements the Negative Binomial (NB) framework to identify significant 
# ribosome stalling sites. It employs a "donut" strategy to calculate local background means by 
# excluding the site of interest, fits a global dispersion trend to account for overdispersion 
# in Ribo-seq count data, and performs FDR-corrected statistical testing to distinguish true 
# translational pauses from stochastic sequencing noise.

# 1. Load Dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("data.table")) install.packages("data.table")
if(!require("MASS")) install.packages("MASS")

suppressPackageStartupMessages({
  library(data.table)
  library(MASS)
  library(stats)
})

#' @details Negative Binomial Fitting and P-value Calculation
fit_nb_stalling_model <- function(final_res, z_window = 200, count_threshold = 10, fc_threshold = 1.5) {
  
  message("Step 1: Calculating 'Donut' Background Means...")
  
  # Ensure data is sorted for rolling window operations
  setorder(final_res, seqnames, pos)
  
  # A. Calculate rolling sums and variances
  # Note: frollsum is used for performance on large genomic data.tables
  final_res[, `:=`(
    temp_sum = frollsum(count, n = z_window, align = "center", fill = NA),
    temp_var = frollapply(count, n = z_window, FUN = var, align = "center", fill = NA)
  )]
  
  # B. Implement Donut Mean Strategy
  # Background Mean = (Window Sum - Current Site Count) / (Window Size - 1)
  # This prevents the signal from masking itself during statistical testing
  final_res[, win_mean := (temp_sum - count) / (z_window - 1)]
  
  # Handle edge cases (division by zero or extremely low background)
  final_res[win_mean <= 0.1, win_mean := 0.1]
  final_res[is.na(win_mean), win_mean := 0.1]
  
  # Recalculate local pause score based on pure background
  final_res[, pause_score := count / win_mean]
  
  # --- Step 2: Fit Global Dispersion Trend ---
  message("Step 2: Fitting Global Dispersion Curve (NLS)...")
  
  # Filter data for robust fitting (avoid low count noise)
  nb_data <- final_res[!is.na(win_mean) & !is.na(temp_var) & win_mean > 1]
  
  # Group by mean bins to stabilize variance estimation
  nb_data[, mean_bin := cut(win_mean, 
                            breaks = quantile(win_mean, probs = seq(0, 1, 0.05)), 
                            include.lowest = TRUE)]
  
  trend_data <- nb_data[, .(
    avg_mean = mean(win_mean, na.rm = TRUE),
    avg_var  = mean(temp_var, na.rm = TRUE)
  ), by = mean_bin]
  
  # Only fit where variance > mean (requirement for Negative Binomial)
  trend_data <- trend_data[avg_var > avg_mean]
  
  # Fit Global Dispersion: Var = Mean + alpha * Mean^2
  fit <- tryCatch({
    nls(avg_var ~ avg_mean + alpha * avg_mean^2, 
        data = trend_data, 
        start = list(alpha = 0.1),
        lower = list(alpha = 0.0001), 
        algorithm = "port")
  }, error = function(e) {
    message("Warning: NLS fit failed, using default dispersion (0.05).")
    return(NULL)
  })
  
  global_alpha <- if (!is.null(fit)) coef(fit)["alpha"] else 0.05
  message(paste("   -> Fitted Global Dispersion (alpha):", round(global_alpha, 4)))
  
  # --- Step 3: Statistical Inference ---
  message("Step 3: Calculating P-values and FDR correction...")
  
  # Define size parameter for NB distribution (1/alpha)
  final_res[, nb_size := 1 / global_alpha]
  
  # Calculate P-values using the upper tail of the NB distribution
  final_res[, p_value := pnbinom(q = count - 1, 
                                 mu = win_mean, 
                                 size = nb_size, 
                                 lower.tail = FALSE)]
  
  # Step 4: Independent Filtering for FDR
  # Focus correction only on sites with sufficient reliability (Candidate Sites)
  candidates_idx <- final_res[count >= count_threshold & pause_score > fc_threshold, which = TRUE]
  
  final_res$padj <- 1.0
  if (length(candidates_idx) > 0) {
    final_res[candidates_idx, padj := p.adjust(p_value, method = "BH")]
  }
  
  return(final_res)
}
