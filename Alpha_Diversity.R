library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(dplyr)

combined_data <- physeq_combined_Algo
metadata_df <- as(sample_data(combined_data), "data.frame")

# ---- Shannon Diversity Analysis ----
alpha_div_shannon_df <- estimate_richness(combined_data, measures = c("Shannon")) %>%
  as.data.frame() %>%
  cbind(metadata_df, .)

alpha_div_shannon_df %>%
  group_by(Project) %>%
  summarise(
    Mean_Shannon = round(mean(Shannon), 2),
    SD_Shannon   = round(sd(Shannon), 2),
    Shannon_Format = paste0(Mean_Shannon, " ± ", SD_Shannon)
  ) %>%
  select(Project, Shannon_Format)

pairwise_comparisons <- list(
  c("Tropical", "Temperate"),
  c("Tropical", "Polar"),
  c("Temperate", "Polar")
)

alpha_div_shannon_plot <- ggplot(alpha_div_shannon_df, aes(x = Geographic_Range, y = Shannon, fill = Geographic_Range)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  labs(x = "Latitude Zones", y = "Shannon Diversity Index", fill = "Latitude Zones") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  scale_fill_viridis_d(option = "D") +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.signif", vjust = 0.5)


anova_shannon <- aov(Shannon ~ Geographic_Range + Location + sample_origin, data = alpha_div_shannon_df)

# ---- Chao1 Diversity Analysis ----
alpha_div_chao1_df <- estimate_richness(combined_data, measures = c("Chao1")) %>%
  as.data.frame() %>%
  cbind(metadata_df, .)

alpha_div_chao1_df %>%
  group_by(Geographic_Range) %>%
  summarise(
    Mean_Chao1 = round(mean(Chao1), 2),
    SD_Chao1   = round(sd(Chao1), 2),
    Chao1_Format = paste0(Mean_Chao1, " ± ", SD_Chao1)
  ) %>%
  select(Geographic_Range, Chao1_Format)

alpha_div_chao1_plot <- ggplot(alpha_div_chao1_df, aes(x = Geographic_Range, y = Chao1, fill = Geographic_Range)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  labs(x = "Latitude Zones", y = "Chao1 Diversity Index", fill = "Latitude Zones") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  scale_fill_viridis_d(option = "D") +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.signif", vjust = 0.5)


anova_chao1 <- aov(Chao1 ~ Geographic_Range + sample_origin + Location, data = alpha_div_chao1_df)

# ---- Number of Genera per Sample ----
otu_table_df <- as.data.frame(otu_table(combined_data))
num_genera <- colSums(otu_table_df > 0)
num_genera_df <- data.frame(SampleID = names(num_genera), Num_Genera = num_genera)
num_genera_df <- merge(metadata_df, num_genera_df, by.x = "row.names", by.y = "SampleID")
colnames(num_genera_df)[1] <- "SampleID"

num_genera_df %>%
  group_by(Geographic_Range) %>%
  summarise(
    Mean_Num_Genera = round(mean(Num_Genera), 2),
    SD_Num_Genera   = round(sd(Num_Genera), 2),
    Num_Genera_Format = paste0(Mean_Num_Genera, " ± ", SD_Num_Genera)
  ) %>%
  select(Geographic_Range, Num_Genera_Format)

num_genera_plot <- ggplot(num_genera_df, aes(x = Geographic_Range, y = Num_Genera, fill = Geographic_Range)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Number of Genera per Sample", x = "Geographic Range", y = "Number of Genera") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_brewer(palette = "Set2") +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.signif")

anova_num_genera <- aov(Num_Genera ~ Geographic_Range + Location + sample_origin, data = num_genera_df)
