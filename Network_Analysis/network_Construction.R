library(SpiecEasi)
library(devtools)
library(igraph)
library(vegan)
library(Matrix)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(future)
library(future.apply)
library(igraph)
library(intergraph)
library(GGally)
library(network)
library(intergraph)
library(RColorBrewer)
library(dplyr)
library(microbiome)


node_gen <- function(merge_data1){
  
  merge_data1.f <- microbiomeutilities::format_to_besthit(merge_data1)
  otu.table <- as(otu_table(merge_data1.f), "matrix")
  env.all <- as(sample_data(merge_data1.f), "data.frame")
  otu.table.all <- t(otu.table)
  row.names(env.all) == row.names(otu.table.all)
  
  if (nrow(otu.table) != nsamples(merge_data1.f)) {
    # If samples are not rows, transpose the OTU table
    otu.table <- t(otu.table)
    
  }
  num_cores <- parallel::detectCores() - 1  # Use one less than the total number of cores to avoid system overload
  run_spiec_easi <- function(otu_table, num_cores) {
    future.apply::future_lapply(1, function(x) {
      spiec.easi(
        otu_table,
        method = 'mb',
        lambda.min.ratio = 1e-2,
        nlambda = 20,
        icov.select.params = list(rep.num = 50, ncores = 1)  # ncores = 1 because future will handle parallelization
      )
    }, future.seed = TRUE)[[1]]  # Set future.seed = TRUE as an argument to future_lapply
  }
  
  # Set up the parallel backend
  spieceasi.net <-  run_spiec_easi(otu.table, num_cores)
  # The result is an object of class "spiec.easi"
  print(spieceasi.net)
  
  # Reset the plan to default
  plan(sequential)
  
  ##Generating the adjacency matrix
  adjacency_matrix <- symBeta(getOptBeta(spieceasi.net))
  colnames(adjacency_matrix) <- rownames(adjacency_matrix) <- colnames(otu.table)
  adjacency_matrix_return <- as.matrix(adjacency_matrix)
  
  return(adjacency_matrix_return)
}




tropical_physeq <- subset_samples(physeq_combined_Algo, 
                                  Geographic_Range %in% c("Tropical"))
temperate_physeq <- subset_samples(physeq_combined_Algo, 
                                   Geographic_Range %in% c("Temperate"))
polar_physeq <- subset_samples(physeq_combined_Algo, 
                               Geographic_Range %in% c("Polar"))
rm(list = setdiff(ls(), c("tropical_physeq", "temperate_physeq", "polar_physeq")))
gc()
complete_data1 <- core(physeq_combined_Algo, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data1 <- prune_samples(non_zero_samples, complete_data1)
complete_data1 <- prune_taxa(taxa_sums(complete_data1) > 0, complete_data1)
adj_complete <- node_gen(complete_data1)

rownames(adj_complete) <- gsub(":.*", "", rownames(adj_complete))
colnames(adj_complete) <- gsub(":.*", "", colnames(adj_complete))

taxa_names <- rownames(adj_complete)

# Subset taxonomy table from phyloseq object
tax_table_df <- as.data.frame(tax_table(complete_data1))
# Make sure rownames match those in adj_complete
common_taxa <- intersect(taxa_names, rownames(tax_table_df))
phyla <- tax_table_df[common_taxa, "Phylum"]
phyla_clean <- phyla[!is.na(phyla)]
phyla_counts <- table(phyla_clean)
phyla_percent <- round(100 * phyla_counts / sum(phyla_counts), 2)

library(ggplot2)

# Create a data frame from phyla_percent
phyla_df <- data.frame(
  Phylum = names(phyla_percent),
  Percentage = as.numeric(phyla_percent)
)

# Remove 'p__' prefix from phylum names
phyla_df$Phylum <- gsub("^p__", "", phyla_df$Phylum)

# Filter to keep only phyla with >1% abundance
phyla_df <- phyla_df[phyla_df$Percentage > 1, ]
# Clean phylum names: remove p__ prefix and anything after first underscore
phyla_df$Phylum <- gsub("_.*", "", gsub("^p__", "", phyla_df$Phylum))
p <- ggplot(phyla_df, aes(x = reorder(Phylum, -Percentage), y = Percentage, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(title = "",
       x = "Phylum",
       y = "Percentage") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size=14),
    axis.text.y = element_text(size=14),
    legend.position = "none"
  ) +
  scale_fill_manual(values = scales::hue_pal()(length(phyla_df$Phylum)))
# Save the plot
ggsave("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/Rendered/phyla_barplot.png", plot = p, width = 8, height = 6, dpi = 800, bg = "transparent")


tropical_physeq1 <- core(tropical_physeq, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(tropical_physeq1) > 0
tropical_physeq1 <- prune_samples(non_zero_samples, tropical_physeq1)
tropical_physeq1 <- prune_taxa(taxa_sums(tropical_physeq1) > 0, tropical_physeq1)
adj_tropical <- node_gen(tropical_physeq1)

temperate_physeq1 <- core(temperate_physeq, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(temperate_physeq1) > 0
temperate_physeq1 <- prune_samples(non_zero_samples, temperate_physeq1)
temperate_physeq1 <- prune_taxa(taxa_sums(temperate_physeq1) > 0, temperate_physeq1)
adj_temperate <- node_gen(temperate_physeq1)

polar_physeq1 <- core(polar_physeq, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(polar_physeq1) > 0
polar_physeq1 <- prune_samples(non_zero_samples, polar_physeq1)
polar_physeq1 <- prune_taxa(taxa_sums(polar_physeq) > 0, polar_physeq)
adj_polar <- node_gen(polar_physeq1)

combined <- netConstruct(data=adj_complete,
                         normMethod = "none", zeroMethod = "none",
                         sparsMethod = "none", dataType = "condDependence",
                         verbose = 3)

adj_tropical <- as.matrix(read.csv("C:/Desktop/PranathiR/adj_tropical.csv", row.names = 1, check.names = FALSE))
adj_temperate <- as.matrix(read.csv("C:/Desktop/PranathiR/adj_temperate.csv", row.names = 1, check.names = FALSE))
adj_polar <- as.matrix(read.csv("C:/Desktop/PranathiR/adj_polar.csv", row.names = 1, check.names = FALSE))

marine_net <- netConstruct(data=adj_marine,
                           normMethod = "none", zeroMethod = "none",
                           sparsMethod = "none", dataType = "condDependence",
                           verbose = 3)

marine_analysis <- netAnalyze(marine_net,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                              hubPar = "degree", hubQuant = 0.95,
                              normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)

sediment_net <- netConstruct(data=adj_sediment,
                             normMethod = "none", zeroMethod = "none",
                             sparsMethod = "none", dataType = "condDependence",
                             verbose = 3)

sediment_analysis <- netAnalyze(sediment_net,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                hubPar = "degree", hubQuant = 0.95,
                                normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)

coral_sponge_net <- netConstruct(data=adj_coral_sponge,
                                 normMethod = "none", zeroMethod = "none",
                                 sparsMethod = "none", dataType = "condDependence",
                                 verbose = 3)

coral_sponge_analysis <- netAnalyze(coral_sponge_net,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                    hubPar = "degree", hubQuant = 0.95,
                                    normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)



