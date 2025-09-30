library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(NetCoMi)
library(microbiome)
library(RColorBrewer)
library(ape)
library(ggrepel)

# === Data Cleaning ===
remove_prefixes <- function(physeq_1) {
  tax_table1 = tax_table(physeq_1)
  taxa_names <- gsub("^d_|^p_|^c_|^o_|^f_|^s_|^_", "", tax_table1)
  tax_names2 <- tax_table(taxa_names)
  otu_table1 <- otu_table(physeq_1)
  metadata <- sam_data(physeq_1)
  phyloseq(otu_table1, tax_names2, metadata)
}

combined_data1 <- remove_prefixes(physeq_combined_Algo)
combined_genus_agglomeration1 <- combined_data1
complete_data1 <- prune_taxa(taxa_sums(combined_genus_agglomeration1) > 0, combined_genus_agglomeration1)
complete_data1 <- core(complete_data1, detection=50, prevalence=0.001, include.lowest=FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data1 <- prune_samples(non_zero_samples, complete_data1)
combined_N <- transform_sample_counts(complete_data1, function(x) x / sum(x))

# === Beta Diversity and PCoA ===
merged_PCoA_N <- ordinate(combined_N, method="PCoA", distance="bray")
bray_dist <- phyloseq::distance(combined_N, method="bray")
meta <- meta(combined_N)

# === PERMDISP (Beta Dispersion) ===
location_betadisp <- betadisper(bray_dist, meta$Geographic_Range)
groups <- location_betadisp$group
betadisp_colors <- brewer.pal(length(unique(groups)), "Set2")
betadisp_test <- permutest(location_betadisp, permutations=999, group=meta$Geographic_Range, pairwise=TRUE)

disp_df <- data.frame(
  Sample = rownames(as.data.frame(location_betadisp$distances)),
  Geographic_Range = location_betadisp$group,
  Distance_to_Centroid = location_betadisp$distances
)

pairwise_comparisons <- list(c("Tropical", "Temperate"), c("Tropical", "Polar"), c("Temperate", "Polar"))

beta_disp_plot <- ggplot(disp_df, aes(x=Geographic_Range, y=Distance_to_Centroid, fill=Geographic_Range)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  labs(x="Latitude Zones", y="Distance to Centroid (Dispersion)", fill="Latitude Zones") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
        axis.text.y=element_text(size=14),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.position="none",
        panel.border=element_rect(color="black", fill=NA, linewidth=0.8)) +
  scale_fill_viridis_d(option="D") +
  stat_compare_means(comparisons=pairwise_comparisons, method="wilcox.test", label="p.signif", vjust=0.5)

ggsave(".../Beta_Dispersion_Geo_Wilcox_colorblind_violin2.png", plot=beta_disp_plot, width=8, height=6, dpi=800)

# === PERMDISP Summary Table ===
group_means <- tapply(location_betadisp$distances, groups, mean, na.rm=TRUE)
group_vars <- tapply(location_betadisp$distances, groups, var, na.rm=TRUE)
perm_results <- betadisp_test
overall_p <- perm_results$tab$`Pr(>F)`[1]
null_expectation <- mean(location_betadisp$distances, na.rm=TRUE)

f_stats <- sapply(1:length(group_means), function(i) group_vars[i] / mean(group_vars[-i], na.rm=TRUE))
p_values <- rep(overall_p, length(group_means))

permdisp_table <- data.frame(
  `Dissimilarity of actual communities` = round(group_means, 3),
  `Dissimilarity of the null expectations` = round(rep(null_expectation, length(group_means)), 3),
  `F` = round(f_stats, 2),
  `p` = ifelse(p_values < 0.001, "< 0.001***", 
               ifelse(p_values < 0.01, paste0(round(p_values, 3), "**"),
                      ifelse(p_values < 0.05, paste0(round(p_values, 3), "*"),
                             round(p_values, 3)))),
  row.names = names(group_means)
)

print(permdisp_table)

# === PERMANOVA ===
adonis_result <- adonis2(bray_dist ~ Geographic_Range*Location+sample_origin, 
                         data=as(sample_data(combined_N), "data.frame"),
                         permutations=999, method="bray", by="terms")
adonis_result

# === Scree Plot ===
princoor <- pcoa(bray_dist)
eigenvalues <- data.frame(PC=1:10, Relative_Eig=princoor$values$Relative_eig[1:10])

scree_plot <- ggplot(eigenvalues, aes(x=factor(PC), y=Relative_Eig)) +
  geom_bar(stat="identity", fill="skyblue", color="black") +
  theme_minimal() +
  labs(x="Principal Coordinate", y="Fraction of Variance Explained") +
  theme(text=element_text(size=14), plot.title=element_text(hjust=0.5))

ggsave(".../scree_plot.png", scree_plot, width=6, height=4, dpi=800)

# === ANOSIM ===
geographi_range.anosim <- anosim(bray_dist, meta$Geographic_Range)
sample_origin_anosim <- anosim(bray_dist, meta$sample_origin)
location_anosim <- anosim(bray_dist, meta$Location)
project_anosim <- anosim(bray_dist, meta$Project)

# === Canonical Correspondence Analysis (CCA) ===
cca_result <- ordinate(combined_N, method="CCA", formula=~ Geographic_Range + Location + sample_origin)
env_scores <- as.data.frame(scores(cca_result, display="bp"))
env_scores$Variable <- gsub("Geographic_Range|sample_origin|Location", "", rownames(env_scores))

arrow_scale <- 2
cca_plot <- plot_ordination(combined_N, cca_result, type="samples", color="Geographic_Range") +
  scale_color_viridis_d(option="D") +
  geom_point(size=1, alpha=0.8) +
  geom_segment(data=env_scores, aes(x=0, y=0, xend=CCA1*arrow_scale*1.5, yend=CCA2*arrow_scale*1.5),
               arrow=arrow(length=unit(0.3, "inches")), color="blue", linewidth=1.0) +
  geom_text_repel(data=env_scores, aes(x=CCA1*arrow_scale*1.6, y=CCA2*arrow_scale*1.6, label=Variable),
                  color="blue", size=3, fontface="bold", max.overlaps=10) +
  labs(color="Latitude Zones", x="CCA1", y="CCA2") +
  theme_bw() +
  theme(axis.text=element_text(color="black"),
        axis.title=element_text(color="black", size=14),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.8))

ggsave(".../cca_plot_colorblind2.png", plot=cca_plot, width=8, height=6, dpi=800)

# === Variation Partitioning Analysis (VPA) ===
otu_data <- as(otu_table(combined_N), "matrix")
otu_data <- t(otu_data)
meta <- data.frame(sample_data(combined_N))

vpa_result <- varpart(bray_dist, ~Location, ~Geographic_Range, ~sample_origin, data=meta)

png(".../vpa_plot.png", width=3600, height=3200, res=300)
plot(vpa_result, digits=2, Xnames=c("Location", "Latitude Zones", "Sample Origin", "Project"))
dev.off()

# === CCA Tables with Envfit ===
meta_clean <- meta[, colSums(is.na(meta)) == 0]
meta_clean$Location_code <- as.numeric(meta_clean$Location_code)
meta_clean$Geographic_Range_code <- as.numeric(meta_clean$Geographic_Range_code)
meta_clean$sample_origin_code <- as.numeric(meta_clean$sample_origin_code)

cca_result_2 <- ordinate(combined_N, method="CCA", formula=~ Geographic_Range_code * Location_code + sample_origin)
env_fit <- envfit(cca_result_2, meta_clean, permutations=999)

cca_results <- data.frame(
  Variable=rownames(env_fit$vectors$arrows),
  CCA1=env_fit$vectors$arrows[,1],
  CCA2=env_fit$vectors$arrows[,2],
  r2=env_fit$vectors$r,
  p_value=env_fit$vectors$pvals
) %>% mutate(Significance=ifelse(p_value < 0.05, "**", ""))