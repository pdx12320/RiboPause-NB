# @title RiboPause-NB: Genomic Context & Feature Extraction
# @description This module extracts localized nucleotide sequences and continuous pause-score tracks 
# surrounding detected stalling sites. It incorporates an automated padding mechanism to ensure 
# uniform feature lengths, facilitating the construction of high-quality input matrices for 
# downstream machine learning and deep learning architectures.

# 1. Load Dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_packages <- c("Rsamtools", "GenomicAlignments", "Biostrings", "data.table")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(Biostrings)
  library(data.table)
})

#' @details Extract padded sequences and score tracks for ML features
extract_context_features <- function(final_res, fa_file, psite_counts, 
                                     up_len = 50, down_len = 50, 
                                     map_clean_to_raw = NULL) {
  
  message("Step 1: Preparing Genomic Tracks...")
  
  # Load FASTA and Index
  fa <- FaFile(fa_file)
  fa_idx <- scanFaIndex(fa)
  fa_lengths <- seqlengths(fa_idx)
  
  # Create a coverage track of pause scores or counts across the transcriptome
  # This serves as the source for 'downstream_scores' and 'upstream_scores'
  score_track <- coverage(psite_counts)
  
  # Handle ID mapping (BAM IDs often lack the version suffix found in FASTA)
  if (is.null(map_clean_to_raw)) {
    raw_ids <- names(fa_lengths)
    clean_ids <- sub("\\..*", "", raw_ids)
    map_clean_to_raw <- setNames(raw_ids, clean_ids)
  }
  
  message("Step 2: Extracting Padded Sequences...")
  
  # Map result IDs to FASTA IDs
  real_fasta_ids <- map_clean_to_raw[final_res$seqnames]
  
  # Define Upstream and Downstream Ranges
  up_gr <- GRanges(real_fasta_ids, IRanges(final_res$pos - up_len, final_res$pos - 1))
  down_gr <- GRanges(real_fasta_ids, IRanges(final_res$pos, final_res$pos + down_len))
  
  # Set sequence lengths and trim to stay within boundaries
  seqlengths(up_gr) <- fa_lengths[seqlevels(up_gr)]
  seqlengths(down_gr) <- fa_lengths[seqlevels(down_gr)]
  up_gr <- trim(up_gr)
  down_gr <- trim(down_gr)
  
  # Fetch sequences
  final_res$upstream_seq <- as.character(getSeq(fa, up_gr))
  final_res$downstream_seq <- as.character(getSeq(fa, down_gr))
  
  message("Step 3: Extracting and Padding Score Tracks...")
  
  # Helper function to extract tracks and enforce fixed lengths via padding
  extract_padded_scores <- function(track, res_dt, target_len, is_upstream = TRUE) {
    
    # Define relative ranges based on position
    if (is_upstream) {
      starts <- res_dt$pos - target_len
      ends <- res_dt$pos - 1
    } else {
      starts <- res_dt$pos
      ends <- res_dt$pos + target_len
    }
    
    # Iterate through sequences to extract views
    all_scores <- vector("list", nrow(res_dt))
    unique_seqs <- unique(res_dt$seqnames)
    
    for (seq in unique_seqs) {
      idx <- which(res_dt$seqnames == seq)
      if (!(seq %in% names(track))) next
      
      # Extract raw values from Rle track
      # Pmax/Pmin ensures we stay within transcript bounds before padding
      chr_len <- length(track[[seq]])
      valid_starts <- pmax(1, starts[idx])
      valid_ends <- pmin(chr_len, ends[idx])
      
      v <- Views(track[[seq]], IRanges(valid_starts, valid_ends))
      val_list <- viewApply(v, as.numeric)
      
      # Apply padding to ensure every entry is exactly target_len
      padded_vals <- sapply(val_list, function(x) {
        if (length(x) < target_len) {
          # Right-pad with zeros for consistency
          x <- c(x, rep(0, target_len - length(x)))
        } else if (length(x) > target_len) {
          x <- x[1:target_len]
        }
        paste(round(x, 3), collapse = ",")
      })
      
      all_scores[idx] <- padded_vals
    }
    return(unlist(all_scores))
  }
  
  final_res$upstream_scores <- extract_padded_scores(score_track, final_res, up_len, TRUE)
  final_res$downstream_scores <- extract_padded_scores(score_track, final_res, down_len, FALSE)
  
  return(final_res)
}
