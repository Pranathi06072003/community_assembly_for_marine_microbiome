# --- Libraries ---
library(iCAMP)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(phyloseq)

# --- Data Preprocessing ---
physeq_clean <- subset_samples(physeq_combined_Algo, sample_origin == "marine metagenome")
physeq_clean <- remove_prefixes(physeq_clean)

physeq_clean <- prune_taxa(taxa_sums(physeq_clean) > 0, physeq_clean)
physeq_clean <- core(physeq_clean, detection = 50, prevalence = 0.01)
physeq_clean <- prune_samples(sample_sums(physeq_clean) > 0, physeq_clean)

# --- Community + Metadata ---
comm <- as.matrix(otu_table(physeq_clean))
if (taxa_are_rows(physeq_clean)) comm <- t(comm)

meta <- data.frame(sample_data(physeq_clean))
treat_vec <- data.frame(meta$Geographic_Range)
row.names(treat_vec) <- sample_names(physeq_clean)

# --- iCAMP Neutral Model ---
neutral_taxa_result <- snm.comm(
  comm, 
  treat = treat_vec, 
  rand = 199,
  alpha = 0.05, 
  two.tail = TRUE,
  output.detail = TRUE
)

# --- Results Formatting ---
df_ratio <- neutral_taxa_result$ratio.summary
df_ratio$treatment.id <- factor(df_ratio$treatment.id,
                                levels = c("Polar", "Temperate", "Tropical"))

custom_colors <- c(
  Polar = "#66c2a5",      # mint green
  Temperate = "#fc8d62",  # orange
  Tropical = "#8da0cb"    # purple-blue
)

# --- Plot ---
p <- ggboxplot(
  df_ratio, 
  x = "treatment.id", 
  y = "Neutral.uw",
  fill = "treatment.id", 
  add = "jitter",
  outlier.shape = NA
) +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "Latitude Zone",
    y = "Proportion of Neutral Taxa"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(p)

# --- Save ---
ggsave("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/Rendered/neutral_taxa_boxplot_marine_prev10.png",
       plot = p, width = 8, height = 6, dpi = 800)
