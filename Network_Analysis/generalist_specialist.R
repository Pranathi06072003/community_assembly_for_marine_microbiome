##========================
## Load required libraries
##========================
library(SpiecEasi)
library(devtools)
library(igraph)
library(vegan)
library(Matrix)
library(reshape2)
library(plyr)
library(future)
library(future.apply)
library(dplyr)
library(gridExtra)
library(grid)
library(intergraph)
library(GGally)
library(network)
library(RColorBrewer)
library(microbiome)
library(NetCoMi)
library(readr)
library(ggplot2)

##========================
## Prepare OTU and Phyloseq data
##========================
complete_data2 <- complete_data1

# Normalize counts
normalized_physeq <- transform_sample_counts(complete_data2, function(x) x / sum(x))

# OTU occupancy
otu_df <- as.data.frame(otu_table(complete_data2))
otu_df$Occupancy <- rowSums(otu_df > 0)

# Define generalists and specialists
threshold_generalist <- 0.50 * ncol(otu_df)
threshold_specialist <- 0.05 * ncol(otu_df)

generalists <- otu_df %>% filter(Occupancy >= threshold_generalist)
generalist_names <- rownames(generalists)

specialists <- otu_df %>% filter(Occupancy <= threshold_specialist)
specialist_taxa <- rownames(specialists)

##========================
## Construct network
##========================
# Remove negative edges
adj_complete[adj_complete < 0] <- 0

network_igraph <- graph_from_adjacency_matrix(adj_complete, mode = "undirected", weighted = TRUE)

# Clean node names
V(network_igraph)$name <- sub(":.*", "", V(network_igraph)$name)

# Compute degrees
degree_list <- degree(network_igraph)
names(degree_list) <- sub(":.*", "", names(degree_list))

generalist_degrees <- degree(network_igraph, v = generalist_names)
specialist_degrees <- degree(network_igraph, v = specialist_taxa)

avg_generalist_degree <- mean(generalist_degrees, na.rm = TRUE)
avg_specialist_degree <- mean(specialist_degrees, na.rm = TRUE)

# Degree distributions
generalist_dist <- table(generalist_degrees) / length(generalist_degrees)
specialist_dist <- table(specialist_degrees) / length(specialist_degrees)

# Plot degree distributions
plot(as.numeric(names(specialist_dist)), as.numeric(specialist_dist), type = "o", pch = 1, col = "red", lwd = 2,
     xlab = "Degree", ylab = "Frequency", main = "Degree Distributions",
     ylim = c(0, max(c(as.numeric(generalist_dist), as.numeric(specialist_dist)))))
points(as.numeric(names(generalist_dist)), as.numeric(generalist_dist), type = "o", pch = 2, col = "blue", lwd = 2)
legend("topright", legend = c("Generalists", "Specialists"), col = c("blue", "red"), pch = c(2,1), lwd = 2)

##========================
## Intra-phyla and Inter-phyla interactions
##========================
taxa_table <- as.data.frame(tax_table(complete_data2))
taxa_to_phylum <- setNames(taxa_table$Phylum, rownames(taxa_table))

edges <- as.data.frame(get.edgelist(network_igraph))
colnames(edges) <- c("Node1", "Node2")

# Filter edges for generalists and specialists
generalist_edges <- edges[edges$Node1 %in% generalist_names | edges$Node2 %in% generalist_names, ]
specialist_edges <- edges[edges$Node1 %in% specialist_taxa | edges$Node2 %in% specialist_taxa, ]

# Add Phylum info
generalist_edges$Phylum1 <- taxa_to_phylum[generalist_edges$Node1]
generalist_edges$Phylum2 <- taxa_to_phylum[generalist_edges$Node2]

specialist_edges$Phylum1 <- taxa_to_phylum[specialist_edges$Node1]
specialist_edges$Phylum2 <- taxa_to_phylum[specialist_edges$Node2]

# Count interactions
intra_phyla_generalists <- sum(generalist_edges$Phylum1 == generalist_edges$Phylum2, na.rm = TRUE)
inter_phyla_generalists <- sum(generalist_edges$Phylum1 != generalist_edges$Phylum2, na.rm = TRUE)

intra_phyla_specialists <- sum(specialist_edges$Phylum1 == specialist_edges$Phylum2, na.rm = TRUE)
inter_phyla_specialists <- sum(specialist_edges$Phylum1 != specialist_edges$Phylum2, na.rm = TRUE)

# Interaction matrix and statistical test
interaction_matrix <- matrix(c(
  intra_phyla_generalists, inter_phyla_generalists,
  intra_phyla_specialists, inter_phyla_specialists
), nrow = 2, byrow = TRUE)
colnames(interaction_matrix) <- c("Intra-phyla", "Inter-phyla")
rownames(interaction_matrix) <- c("Generalists", "Specialists")

test_result <- if (any(interaction_matrix < 5)) fisher.test(interaction_matrix) else chisq.test(interaction_matrix)
print(test_result)

# Bar plot of interactions
interaction_df <- data.frame(
  Group = rep(c("Generalists", "Specialists"), each = 2),
  Interaction_Type = rep(c("Intra-phyla", "Inter-phyla"), 2),
  Count = c(intra_phyla_generalists, inter_phyla_generalists, intra_phyla_specialists, inter_phyla_specialists)
)
interaction_plot <- ggplot(interaction_df, aes(x = Group, y = Count, fill = Interaction_Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.3) +
  labs(x = "Taxa Group", y = "Number of Interactions", fill = "Interaction Type") +
  theme_minimal()
ggsave("../interaction_plot_generalists.png",
       plot = interaction_plot, width = 6, height = 4, dpi = 800)

##========================
## Phylum distribution significance
##========================
phyla_counts_g <- generalist_taxa %>% dplyr::count(Phylum)
phyla_counts_s <- specialist_taxa2 %>% dplyr::count(Phylum)

colnames(phyla_counts_g) <- c("Phylum", "n_G")
colnames(phyla_counts_s) <- c("Phylum", "n_S")

phyla_combined <- merge(phyla_counts_g, phyla_counts_s, by = "Phylum", all = TRUE)
phyla_combined[is.na(phyla_combined)] <- 0
phyla_combined$n_G <- as.numeric(phyla_combined$n_G)
phyla_combined$n_S <- as.numeric(phyla_combined$n_S)

# Fisher test per phylum
phyla_combined$p_value <- apply(phyla_combined, 1, function(row) {
  contingency <- matrix(c(as.numeric(row["n_G"]), sum(phyla_counts_g$n_G) - as.numeric(row["n_G"]),
                          as.numeric(row["n_S"]), sum(phyla_counts_s$n_S) - as.numeric(row["n_S"])),
                        nrow = 2, byrow = TRUE)
  fisher.test(contingency)$p.value
})
phyla_combined$log_p <- -log10(phyla_combined$p_value)
phyla_combined$Significance <- ifelse(phyla_combined$p_value < 0.05, "Significant", "Not Significant")
phyla_combined$Phylum <- sub("p__*", "", phyla_combined$Phylum)

# Plot significance
ggplot(phyla_combined %>% filter(log_p > 0 & Phylum != "SAR324"), 
       aes(x = reorder(Phylum, log_p), y = log_p, fill = Significance)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Significant" = "darkblue", "Not Significant" = "grey")) +
  labs(x = "Phylum", y = "-log(p)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
ggsave("../phylum_significance_plot.jpg",
       width = 8, height = 6, dpi = 800)

##========================
## Zi-Pi analysis (module connectivity)
##========================
igraph.group1 <- graph_from_adjacency_matrix(adj_complete, mode = "undirected", weighted = TRUE, diag = FALSE)
igraph.group1.weight <- E(igraph.group1)$weight
E.color.group1 <- ifelse(igraph.group1.weight > 0, "pink", ifelse(igraph.group1.weight < 0, "blue", "grey"))
E(igraph.group1)$color <- E.color.group1
E(igraph.group1)$width <- abs(igraph.group1.weight)

fc.group1 <- cluster_fast_greedy(igraph.group1)
modularity.group1 <- modularity(fc.group1)
Membership <- membership(fc.group1)

# Calculate Zi and Pi
Nnodes <- vcount(igraph.group1)
Z <- P <- numeric(Nnodes)
for (i in 1:Nnodes) {
  L <- Membership == Membership[i]
  neighbs <- neighbors(igraph.group1, i)
  Kis <- sum(L[neighbs])
  SUM <- SUMsq <- SUMP <- 0
  Miv <- which(L)
  for (j in Miv) {
    neighbsj <- neighbors(igraph.group1, j)
    Kjs <- sum(L[neighbsj])
    SUM <- SUM + Kjs
    SUMsq <- SUMsq + Kjs^2
  }
  Z[i] <- ifelse((Kis - SUM / sum(L)) == 0, 0, (Kis - SUM / sum(L)) / sqrt(SUMsq / sum(L) - (SUM / sum(L))^2))
  
  for (k in unique(Membership)) {
    Lp <- Membership == k
    Kisp <- sum(Lp[neighbs])
    SUMP <- SUMP + (Kisp / degree(igraph.group1, i))^2
  }
  P[i] <- 1 - SUMP
}

# Zi-Pi plot
png("../inter_intra_module_connectivity.png",
    width = 4800, height = 3600, res = 800)
par(mar = c(5, 5, 4, 5), mgp = c(1.8, 0.5, 0))
plot(P, Z, xlim = c(0, 1), ylim = c(-4, 4), xlab = "Among-module connectivity P", ylab = "Within-module connectivity Z",
     col = ifelse(V(igraph.group1)$name %in% generalist_names, "blue", ifelse(V(igraph.group1)$name %in% specialist_taxa, "red", "grey")),
     pch = 19, cex = 0.95)
abline(v = 0.62, h = 2.5, col = "black")
text(0.15, 4, "Module hubs"); text(0.8, 4, "Network hubs")
text(0.15, -4, "Peripherals"); text(0.8, -4, "Connectors")
dev.off()

##========================
## Natural connectivity analysis
##========================
calc_natural_connectivity <- function(graph) {
  adj_matrix <- as.matrix(as_adjacency_matrix(graph, attr = "weight", sparse = FALSE))
  mean(exp(eigen(adj_matrix)$values))
}

removal_percentages <- seq(5, 60, by = 5)
natural_connectivities <- numeric(length(removal_percentages))
igraph_copy <- igraph.group1

for (idx in seq_along(removal_percentages)) {
  percentage <- removal_percentages[idx]
  num_nodes_to_remove <- floor(vcount(igraph_copy) * (percentage / 100))
  
  for (i in 1:num_nodes_to_remove) {
    if (vcount(igraph_copy) == 0) break
    highest_degree_node <- V(igraph_copy)[which.max(degree(igraph_copy))]
    igraph_copy <- delete_vertices(igraph_copy, highest_degree_node)
  }
  natural_connectivities[idx] <- calc_natural_connectivity(igraph_copy)
}

# Plot natural connectivity
plot(removal_percentages, natural_connectivities, type = "b", pch = 19, col = "blue",
     xlab = "Proportion of removed nodes (%)", ylab = "Natural Connectivity",
     main = "Impact of Node Removal on Natural Connectivity")

