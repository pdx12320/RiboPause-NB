# RiboPause-NB

**RiboPause-NB** is a robust R-based pipeline designed for high-resolution ribosome stalling detection using Negative Binomial (NB) statistical modeling and optimized Z-scores. It is specifically engineered to generate ML-ready datasets by extracting precise codon-level sequences and localized pause score tracks with automated padding.

---

## 🚀 Key Features

* **Statistical Precision**: Utilizes Negative Binomial modeling to handle the overdispersion common in Ribo-seq count data.
* **Donut Background Strategy**: Implements a "donut" mean strategy for local background calculation, ensuring the site being tested does not mask its own signal.
* **Machine Learning Integration**: Automatically extracts upstream/downstream sequences and continuous pause-score tracks, using auto-padding to ensure uniform input lengths for downstream ML models.
* **Expression-Aware Filtering**: Built-in logic to filter low-expression genes (e.g., keeping the top 50% by mean depth) to eliminate false positives in low-coverage regions.
* **Version-Agnostic Mapping**: Automatically handles Ensembl ID version mismatches (e.g., cleaning `.1` suffixes) between BAM files and FASTA references.
* **In-frame Annotation**: Filters pause sites to ensure they are in-frame and maps them to specific codons using GTF annotation.

---

## 🛠 Methodology

The pipeline follows a rigorous statistical workflow:

1.  **P-site Mapping**: Reads BAM alignments and applies length-specific offsets to identify precise P-site coordinates.
2.  **Dispersion Fitting**: Fits a global dispersion trend using Non-linear Least Squares (NLS) to model the variance-mean relationship:
    $$Var = \mu + \alpha \cdot \mu^2$$
3.  **Local Z-score Calculation**: Computes a modified Z-score to suppress noise in low-expression regions:
    $$Z = \frac{x - \mu}{\sigma + 1}$$
4.  **P-value & FDR**: Calculates P-values using the Negative Binomial distribution and applies Benjamini-Hochberg (BH) correction only to high-confidence candidate sites.
5.  **Feature Extraction**: Generates serialized context scores for the upstream and downstream regions (e.g., 50nt) with forced padding for consistency.

---

## 📦 Installation

Ensure you have R (>= 4.0) and the following Bioconductor packages installed:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsamtools", "GenomicAlignments", "Biostrings", "GenomicRanges"))
install.packages(c("data.table", "optparse", "MASS"))
```
---

## 💻 Quick Start

You can run the main analysis via the command line:

```
Rscript main.R \
  --bam ./ribo-data/sorted.bam \
  --ref ./ribo-data/Mus_musculus.GRCm39.cdna.all.fa \
  --annot ./ribo-data/Mus_musculus.GRCm39.104.gtf \
  --lengths "28,29,30" \
  --offset "12,12,12" \
  --foldchange 10 \
  --output result.csv
```

## 📊 Output Format
The primary output is a CSV file containing:

* **seqnames/pos**: Transcript ID and P-site coordinate.

* **pause_score**: Fold-change relative to local background.

* **z_score**: The optimized local Z-score.

* **upstream_seq / downstream_seq**: Nucleotide sequences for context.

* **upstream_scores / downstream_scores**: Comma-separated pause score tracks for ML features.

* **codon**: The specific codon at the pause site (if annotation is provided).

## 👤 Author
Yang Dengxiang

Bioinformatics Undergraduate, Southern University of Science and Technology (SUSTech)
