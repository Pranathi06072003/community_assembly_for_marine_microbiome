# ===============================================================
# Required Packages
# ===============================================================
library(phyloseq)
library(vegan)
library(ggplot2)

# ===============================================================
# Function: Remove Taxonomy Prefixes
# ===============================================================
remove_prefixes <- function(physeq_1) {
  tax_table1 = tax_table(physeq_1)
  taxa_names <- gsub("^d_", "", tax_table1)
  taxa_names <- gsub("^p_", "", taxa_names)
  taxa_names <- gsub("^c_", "", taxa_names)
  taxa_names <- gsub("^o_", "", taxa_names)
  taxa_names <- gsub("^f_", "", taxa_names)
  taxa_names <- gsub("^s_", "", taxa_names)
  taxa_names <- gsub("^_", "", taxa_names)
  tax_names2 <- tax_table(taxa_names)
  otu_table1 <- otu_table(physeq_1)
  metadata <- sam_data(physeq_1)
  physeq2 <- phyloseq(otu_table1, tax_names2, metadata)
  return(physeq2)
}

# ===============================================================
# Preprocessing
# ===============================================================
combined_data1 <- remove_prefixes(physeq_combined_Algo)
combined_genus_agglomeration1 <- combined_data1
complete_data1 <- prune_taxa(taxa_sums(combined_genus_agglomeration1) > 0, combined_genus_agglomeration1)
complete_data1 <- core(complete_data1, detection = 50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data1 <- prune_samples(non_zero_samples, complete_data1)
combined_N <- transform_sample_counts(complete_data1, function(x) x / sum(x))

# ===============================================================
# Canonical Correspondence Analysis (CCA)
# ===============================================================
cca_result <- ordinate(combined_N, method = "CCA", 
                       formula = ~ sample_origin + Project + longitude_deg + latitude_deg)

plot_ordination(combined_N, cca_result, color = "sample_origin", 
                title = "Canonical Correlation Analysis After Normalisation")

# ===============================================================
# NMDS with Environmental Fitting (envfit)
# ===============================================================
ord_scores <- ordinate(combined_N, method = "NMDS", distance = "euclidean")
env_data <- metadata[, c("Project", "sample_origin", "longitude_deg", "latitude_deg", "temperature_deg_c")]
envfit_result <- envfit(ord_scores, env_data, permutations = 999, na.rm = TRUE)

env_vectors <- as.data.frame(envfit_result$vectors$arrows)
env_vectors$label <- rownames(env_vectors)
env_vectors$label <- gsub("longitude_deg", "Longitude", env_vectors$label)
env_vectors$label <- gsub("latitude_deg", "Latitude", env_vectors$label)

plot_ordination(combined_N, ord_scores, color = "sample_origin") + 
  geom_point() + 
  geom_segment(data = env_vectors, 
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), 
               color = "blue", size = 1) + 
  ggtitle("NMDS Ordination with Environmental Variables") +
  geom_text(data = env_vectors, 
            aes(x = NMDS1, y = NMDS2, label = label), 
            color = "red", size = 4, check_overlap = TRUE, 
            nudge_x = 0.2, nudge_y = 0.2)

# ===============================================================
# Variation Partitioning Analysis (VPA)
# ===============================================================
physeq <- combined_N
otu_data <- otu_table(physeq)
metadata <- sample_data(physeq)

env1 <- as.data.frame(metadata[, "temperature_deg_c", drop = FALSE])
env2 <- as.data.frame(metadata[, "sample_origin", drop = FALSE])
env3 <- as.data.frame(metadata[, "latitude_deg", drop = FALSE])
env4 <- as.data.frame(metadata[, "longitude_deg", drop = FALSE])

env1_values <- as.data.frame(sample_data(env1))
env2_values <- as.data.frame(sample_data(env2))
env3_values <- as.data.frame(sample_data(env3))
env4_values <- as.data.frame(sample_data(env4))

otu_data <- as.matrix(otu_data)
otu_data <- t(otu_data)
transformed_data <- decostand(otu_data, method = "hellinger")

vpa_result <- varpart(transformed_data, ~ env1_values + env2_values + env3_values + env4_values)
