library(spaa)
library(ggpubr)
library(ggplot2)

all_data <- physeq_combined_Algo

tropical_physeq <- subset_samples(all_data, Geographic_Range %in% c("Tropical"))
temperate_physeq <- subset_samples(all_data, Geographic_Range %in% c("Temperate"))
polar_physeq <- subset_samples(all_data, Geographic_Range %in% c("Polar"))

process_physeq <- function(physeq) {
  physeq_core <- core(physeq, detection = 50, prevalence = 0.01, include.lowest = FALSE)
  non_zero_samples <- sample_sums(physeq_core) > 0
  physeq_core <- prune_samples(non_zero_samples, physeq_core)
  prune_taxa(taxa_sums(physeq_core) > 0, physeq_core)
}

tropical_physeq1 <- process_physeq(tropical_physeq)
temperate_physeq1 <- process_physeq(temperate_physeq)
polar_physeq1 <- process_physeq(polar_physeq)

otu_mat_tropical <- t(as.matrix(otu_table(tropical_physeq1)))
otu_mat_temperate <- t(as.matrix(otu_table(temperate_physeq1)))
otu_mat_polar <- t(as.matrix(otu_table(polar_physeq1)))

result_tropical <- niche.width(otu_mat_tropical, method = "levins")
result_temperate <- niche.width(otu_mat_temperate, method = "levins")
result_polar <- niche.width(otu_mat_polar, method = "levins")

df_tropical <- data.frame(NicheWidth = as.numeric(result_tropical), Climate = "Tropical")
df_temperate <- data.frame(NicheWidth = as.numeric(result_temperate), Climate = "Temperate")
df_polar <- data.frame(NicheWidth = as.numeric(result_polar), Climate = "Polar")

combined_df <- rbind(df_tropical, df_temperate, df_polar)

pairwise_comparisons <- list(
  c("Tropical", "Temperate"),
  c("Tropical", "Polar"),
  c("Temperate", "Polar")
)

p <- ggplot(combined_df, aes(x = Climate, y = NicheWidth, fill = Climate)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0, 350)) +
  labs(x = "Latitude Zone", y = "Niche Width") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.signif")

ggsave(
  "../niche_width_boxplot.png",
  plot = p, width = 8, height = 6, dpi = 800
)

