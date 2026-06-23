###############################################################
# Batch Effect Correction using MMUPHin
# Dataset: Combined Marine Microbiome Dataset
# Method:
#   1. Filter taxa for network analysis
#   2. Normalize abundances
#   3. Remove Project-specific batch effects using MMUPHin
#   4. Compute Bray-Curtis distances
#   5. Evaluate remaining sources of variation using PERMANOVA
###############################################################

# ============================================================
# Load libraries
# ============================================================

library(MMUPHin)
library(phyloseq)
library(vegan)
library(dplyr)

# ============================================================
# 1. Create filtered dataset for downstream analyses
# ============================================================

# Retain taxa present at >=50% abundance in at least 1% of samples
complete_data <- core(
  physeq_combined_Algo,
  detection = 50,
  prevalence = 0.01,
  include.lowest = FALSE
)

# Remove empty samples and taxa
complete_data <- prune_samples(sample_sums(complete_data) > 0, complete_data)
complete_data <- prune_taxa(taxa_sums(complete_data) > 0, complete_data)

# Relative abundance normalization
filtered_physeq <- transform_sample_counts(
  complete_data,
  function(x) x / sum(x)
)

# ============================================================
# 2. Function for MMUPHin batch correction
# ============================================================

run_batch_correction <- function(physeq_object) {
  
  otu_mat <- as(otu_table(physeq_object), "matrix")
  
  metadata <- as(sample_data(physeq_object), "data.frame")
  metadata <- metadata[colnames(otu_mat), , drop = FALSE]
  
  # Convert metadata variables to factors
  metadata$Project_ID <- as.factor(metadata$Project_ID)
  metadata$sample_origin_code <- as.factor(metadata$sample_origin_code)
  metadata$Geographic_Range_code <- as.factor(metadata$Geographic_Range_code)
  metadata$Location_code <- as.factor(metadata$Location_code)
  
  # MMUPHin batch correction
  fit <- adjust_batch(
    feature_abd = otu_mat,
    batch = "Project_ID",
    covariates = NULL,
    data = metadata,
    control = list(verbose = FALSE)
  )
  
  phyloseq(
    otu_table(fit$feature_abd_adj, taxa_are_rows = TRUE),
    tax_table(physeq_object),
    sample_data(physeq_object)
  )
}

# ============================================================
# 3. Batch correction
# ============================================================

# Filtered dataset (used for downstream analyses)
combined_N_NCV <- run_batch_correction(filtered_physeq)

# Full dataset (used elsewhere in the project)
physeq_combined_Algo <- run_batch_correction(physeq_combined_Algo)

# ============================================================
# 4. Bray-Curtis distance matrix
# ============================================================

bray_dist_new_NCV <- phyloseq::distance(
  combined_N_NCV,
  method = "bray"
)

# Uncomment to save
# saveRDS(
#   bray_dist_new_NCV,
#   "Objects_Markdown/bray_dist_new_NCV.rds"
# )

# Uncomment to load precomputed distances
# bray_dist_new_NCV <- readRDS(
#   "Objects_Markdown/bray_dist_new_NCV.rds"
# )

# ============================================================
# 5. PERMANOVA
# ============================================================

sample_data(combined_N_NCV)$target_subfragment <-
  as.factor(sample_data(combined_N_NCV)$target_subfragment)

adonis_result_new_NCV <- adonis2(
  bray_dist_new_NCV ~
    sample_origin +
    Geographic_Range +
    Project +
    Location +
    target_subfragment,
  data = as(sample_data(combined_N_NCV), "data.frame"),
  method = "bray",
  permutations = 999,
  by = "margin",
  parallel = 4
)

# Uncomment to save
# saveRDS(
#   adonis_result_new_NCV,
#   "Objects_Markdown/adonis_result_new_NCV.rds"
# )

# Uncomment to load
# adonis_result_new_NCV <- readRDS(
#   "Objects_Markdown/adonis_result_new_NCV.rds"
# )

# ============================================================
# 6. Results
# ============================================================

adonis_result_new_NCV