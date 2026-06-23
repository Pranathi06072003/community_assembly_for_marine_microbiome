############################################################
## Subsample Tropical Samples
############################################################

set.seed(42)

# Number of samples to retain from each sample origin
n_samples <- 51

# Function to randomly subsample a sample origin
subsample_origin <- function(physeq, origin, zone = "Tropical", n = 51) {
  
  samples <- sample_names(
    subset_samples(
      physeq,
      sample_origin == origin &
        Latitude_Zones == zone
    )
  )
  
  phy <- prune_samples(sample(samples, n), physeq)
  
  phy <- prune_taxa(taxa_sums(phy) > 0, phy)
  phy <- prune_samples(sample_sums(phy) > 0, phy)
  
  return(phy)
}

physeq_sediment <- subsample_origin(
  physeq_combined_Algo,
  "marine sediment metagenome",
  n = n_samples
)

physeq_marine <- subsample_origin(
  physeq_combined_Algo,
  "marine metagenome",
  n = n_samples
)

physeq_coral <- subsample_origin(
  physeq_combined_Algo,
  "coral metagenome",
  n = n_samples
)

physeq_sponge <- subsample_origin(
  physeq_combined_Algo,
  "sponge metagenome",
  n = n_samples
)

############################################################
## Batch Correction (Marine Samples)
############################################################

batch_correct <- function(physeq){
  
  otu_mat <- as(otu_table(physeq), "matrix")
  
  metadata <- as(sample_data(physeq), "data.frame")
  metadata <- metadata[colnames(otu_mat), , drop = FALSE]
  metadata$Project_ID <- factor(metadata$Project_ID)
  
  fit <- adjust_batch(
    feature_abd = otu_mat,
    batch = "Project_ID",
    covariates = NULL,
    data = metadata,
    control = list(verbose = FALSE)
  )
  
  otu_new <- otu_table(
    fit$feature_abd_adj,
    taxa_are_rows = TRUE
  )
  
  phyloseq(
    otu_new,
    tax_table(physeq),
    sample_data(physeq)
  )
}

physeq_marine <- batch_correct(physeq_marine)

############################################################
## Merge Subsampled Datasets
############################################################

physeq_subsampled <- merge_phyloseq(
  physeq_sediment,
  physeq_marine,
  physeq_coral,
  physeq_sponge
)

physeq_subsampled <- prune_taxa(
  taxa_sums(physeq_subsampled) > 0,
  physeq_subsampled
)

############################################################
## PERMANOVA
############################################################

bray_dist_SS <- phyloseq::distance(
  physeq_subsampled,
  method = "bray"
)

adonis_result_SS <- adonis2(
  bray_dist_SS ~ sample_origin +
    Geographic_Range +
    Project +
    Location,
  data = as(sample_data(physeq_subsampled), "data.frame"),
  permutations = 999,
  method = "bray",
  by = "margin",
  parallel = 4
)

adonis_result_SS

############################################################
## Alpha Diversity
############################################################

combined_data <- physeq_subsampled

alpha_div_df <- estimate_richness(
  combined_data,
  measures = "Shannon"
) %>%
  cbind(as(sample_data(combined_data), "data.frame"), .)

############################################################
## Summary Statistics
############################################################

alpha_div_df %>%
  group_by(sample_origin) %>%
  summarise(
    Mean_Shannon = round(mean(Shannon), 2),
    SD_Shannon = round(sd(Shannon), 2),
    Shannon = paste0(Mean_Shannon, " ± ", SD_Shannon),
    .groups = "drop"
  ) %>%
  select(sample_origin, Shannon)

############################################################
## Pairwise Comparisons
############################################################

pairwise_comparisons <- list(
  c("marine metagenome", "marine sediment metagenome"),
  c("marine metagenome", "coral metagenome"),
  c("marine metagenome", "sponge metagenome"),
  c("coral metagenome", "sponge metagenome"),
  c("coral metagenome", "marine sediment metagenome"),
  c("sponge metagenome", "marine sediment metagenome")
)

############################################################
## Shannon Diversity Plot
############################################################

alpha_div_plot <- ggplot(
  alpha_div_df,
  aes(sample_origin, Shannon, fill = sample_origin)
) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    comparisons = pairwise_comparisons,
    method = "wilcox.test",
    label = "p.signif"
  ) +
  scale_fill_viridis_d(option = "D") +
  labs(
    x = "Sample Origin",
    y = "Shannon Diversity Index"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none",
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

alpha_div_plot

############################################################
## Save Figure
############################################################

ggsave(
  "results/Alpha_Diversity_Shannon.png",
  alpha_div_plot,
  width = 8,
  height = 6,
  dpi = 800
)

############################################################
## ANOVA
############################################################

anova_shannon <- aov(
  Shannon ~ Geographic_Range +
    Location +
    sample_origin,
  data = alpha_div_df
)

summary(anova_shannon)

############################################################
## Tukey HSD
############################################################

tukey_shannon <- TukeyHSD(anova_shannon)

print(tukey_shannon)

############################################################
## Save Results
############################################################

write.csv(
  data.frame(summary(anova_shannon)[[1]]),
  "results/Shannon_ANOVA.csv"
)

write.csv(
  as.data.frame(tukey_shannon$Geographic_Range),
  "results/Shannon_TukeyHSD.csv"
)

############################################################
## Chao1 Diversity
############################################################

alpha_div_df <- estimate_richness(
  combined_data,
  measures = "Chao1"
) %>%
  cbind(as(sample_data(combined_data), "data.frame"), .)

############################################################
## Chao1 Diversity
############################################################

alpha_div_df <- estimate_richness(
  combined_data,
  measures = "Chao1"
) %>%
  cbind(as(sample_data(combined_data), "data.frame"), .)

############################################################
## Summary Statistics
############################################################

alpha_div_df %>%
  group_by(sample_origin) %>%
  summarise(
    Mean_Chao1 = round(mean(Chao1, na.rm = TRUE), 2),
    SD_Chao1 = round(sd(Chao1, na.rm = TRUE), 2),
    Chao1 = paste0(Mean_Chao1, " ± ", SD_Chao1),
    .groups = "drop"
  ) %>%
  select(sample_origin, Chao1)

############################################################
## Pairwise Comparisons
############################################################

pairwise_comparisons <- list(
  c("marine metagenome", "marine sediment metagenome"),
  c("marine metagenome", "coral metagenome"),
  c("marine metagenome", "sponge metagenome"),
  c("coral metagenome", "sponge metagenome"),
  c("coral metagenome", "marine sediment metagenome"),
  c("sponge metagenome", "marine sediment metagenome")
)

############################################################
## Chao1 Boxplot
############################################################

alpha_div_plot <- ggplot(
  alpha_div_df,
  aes(sample_origin, Chao1, fill = sample_origin)
) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    comparisons = pairwise_comparisons,
    method = "wilcox.test",
    label = "p.signif"
  ) +
  scale_fill_viridis_d(option = "D") +
  labs(
    x = "Sample Origin",
    y = "Chao1 Diversity Index"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none",
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

alpha_div_plot

############################################################
## Save Figure
############################################################

ggsave(
  "results/Alpha_Diversity_Chao1.png",
  alpha_div_plot,
  width = 8,
  height = 6,
  dpi = 800
)

############################################################
## ANOVA
############################################################

anova_chao1 <- aov(
  Chao1 ~ Geographic_Range +
    sample_origin +
    Location,
  data = alpha_div_df
)

summary(anova_chao1)

############################################################
## Tukey HSD
############################################################

tukey_chao1 <- TukeyHSD(anova_chao1)

print(tukey_chao1)

############################################################
## Save Statistical Results
############################################################

write.csv(
  data.frame(summary(anova_chao1)[[1]]),
  "results/Chao1_ANOVA.csv"
)

write.csv(
  as.data.frame(tukey_chao1$Geographic_Range),
  "results/Chao1_TukeyHSD.csv"
)

############################################################
## Beta Dispersion Analysis
############################################################

# Remove taxonomic prefixes
remove_prefixes <- function(physeq) {
  
  taxonomy <- tax_table(physeq)
  
  taxonomy <- gsub("^[dpcofs]_", "", taxonomy)
  taxonomy <- gsub("^_", "", taxonomy)
  
  phyloseq(
    otu_table(physeq),
    tax_table(taxonomy),
    sample_data(physeq)
  )
}

############################################################
## Core Microbiome Filtering
############################################################

combined_core <- core(
  physeq_subsampled,
  detection = 50,
  prevalence = 0.01,
  include.lowest = FALSE
)

combined_core <- prune_samples(
  sample_sums(combined_core) > 0,
  combined_core
)

combined_core <- prune_taxa(
  taxa_sums(combined_core) > 0,
  combined_core
)

combined_core <- transform_sample_counts(
  combined_core,
  function(x) x / sum(x)
)

############################################################
## Bray-Curtis Distance
############################################################

bray_dist <- phyloseq::distance(
  combined_core,
  method = "bray"
)

metadata <- as(sample_data(combined_core), "data.frame")

saveRDS(
  bray_dist,
  "results/bray_distance_sample_origin.rds"
)

############################################################
## Beta Dispersion
############################################################

location_betadisp <- betadisper(
  bray_dist,
  metadata$sample_origin
)

plot(
  location_betadisp,
  hull = FALSE,
  ellipse = TRUE
)

saveRDS(
  location_betadisp,
  "results/betadisp_sample_origin.rds"
)

############################################################
## Permutation Test
############################################################

betadisp_test <- permutest(
  location_betadisp,
  permutations = 999,
  pairwise = TRUE
)

############################################################
## Distance to Centroid
############################################################

disp_df <- data.frame(
  Sample = names(location_betadisp$distances),
  sample_origin = location_betadisp$group,
  Distance_to_Centroid = location_betadisp$distances
)

############################################################
## Beta Dispersion Summary
############################################################

disp_summary <- disp_df %>%
  group_by(sample_origin) %>%
  summarise(
    Mean_Distance = mean(Distance_to_Centroid, na.rm = TRUE),
    SD_Distance = sd(Distance_to_Centroid, na.rm = TRUE),
    Sample_Count = n(),
    .groups = "drop"
  )

disp_summary

############################################################
## Pairwise Comparisons
############################################################

pairwise_comparisons <- list(
  c("marine metagenome", "marine sediment metagenome"),
  c("marine metagenome", "coral metagenome"),
  c("marine metagenome", "sponge metagenome"),
  c("coral metagenome", "sponge metagenome"),
  c("coral metagenome", "marine sediment metagenome"),
  c("sponge metagenome", "marine sediment metagenome")
)

############################################################
## Beta Dispersion Plot
############################################################

beta_disp_plot <- ggplot(
  disp_df,
  aes(sample_origin, Distance_to_Centroid, fill = sample_origin)
) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    comparisons = pairwise_comparisons,
    method = "wilcox.test",
    label = "p.signif"
  ) +
  scale_fill_viridis_d(option = "D") +
  labs(
    x = "Sample Origin",
    y = "Distance to Centroid"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none",
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

beta_disp_plot

############################################################
## Save Figure
############################################################

ggsave(
  "results/Beta_Dispersion_SampleOrigin.png",
  beta_disp_plot,
  width = 8,
  height = 6,
  dpi = 800
)

############################################################
## Network Analysis
############################################################

rm(list = setdiff(
  ls(),
  c(
    "physeq_coral",
    "physeq_sponge",
    "physeq_sediment",
    "physeq_marine",
    "node_gen",
    "adjust_batch",
    "adonis_result_SS"
  )
))

############################################################
## Helper Function: Core Filtering
############################################################

prepare_phyloseq <- function(physeq) {
  
  physeq <- core(
    physeq,
    detection = 50,
    prevalence = 0.01,
    include.lowest = FALSE
  )
  
  physeq <- prune_samples(
    sample_sums(physeq) > 0,
    physeq
  )
  
  prune_taxa(
    taxa_sums(physeq) > 0,
    physeq
  )
}

############################################################
## Coral Network
############################################################

# Coral samples originate from a single project,
# therefore batch correction was not required.

physeq_coral <- prepare_phyloseq(physeq_coral)

adj_coral <- node_gen(
  physeq_coral,
  "results/network_files/nodes_coral.csv",
  "results/network_files/edges_coral.csv",
  "results/network_files/metrics_coral.csv",
  "results/network_files/degree_coral.png"
)

############################################################
## Sponge Network
############################################################

physeq_sponge <- prepare_phyloseq(physeq_sponge)

adj_sponge <- node_gen(
  physeq_sponge,
  "results/network_files/nodes_sponge.csv",
  "results/network_files/edges_sponge.csv",
  "results/network_files/metrics_sponge.csv",
  "results/network_files/degree_sponge.png"
)

############################################################
## Sediment Network
############################################################

physeq_sediment <- prepare_phyloseq(physeq_sediment)

adj_sediment <- node_gen(
  physeq_sediment,
  "results/network_files/nodes_sediment.csv",
  "results/network_files/edges_sediment.csv",
  "results/network_files/metrics_sediment.csv",
  "results/network_files/degree_sediment.png"
)

############################################################
## Marine Network (Batch Corrected)
############################################################

physeq_marine <- prepare_phyloseq(physeq_marine)

otu_mat <- as(otu_table(physeq_marine), "matrix")

metadata <- as(sample_data(physeq_marine), "data.frame")
metadata <- metadata[colnames(otu_mat), , drop = FALSE]
metadata$Project_ID <- factor(metadata$Project_ID)

fit_adjust_batch <- adjust_batch(
  feature_abd = otu_mat,
  batch = "Project_ID",
  covariates = NULL,
  data = metadata,
  control = list(verbose = FALSE)
)

physeq_marine <- phyloseq(
  otu_table(
    fit_adjust_batch$feature_abd_adj,
    taxa_are_rows = TRUE
  ),
  tax_table(physeq_marine),
  sample_data(physeq_marine)
)

adj_marine <- node_gen(
  physeq_marine,
  "results/network_files/nodes_marine.csv",
  "results/network_files/edges_marine.csv",
  "results/network_files/metrics_marine.csv",
  "results/network_files/degree_marine.png"
)

############################################################
## Save Adjacency Matrices
############################################################

write.csv(adj_coral, "results/adjacency_coral.csv")

write.csv(adj_sponge, "results/adjacency_sponge.csv")

write.csv(adj_sediment, "results/adjacency_sediment.csv")

write.csv(adj_marine, "results/adjacency_marine.csv")

##ZiPi for network analyses: 
library(igraph) 
library(psych)

#"group1" is the OTU table, has been filtered to include >1% and presence frequency >2 samples.
#occor.group1 = corr.test(t(group1),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
#occor.r.group1= occor.group1$r 
#occor.p.group1= occor.group1$p 
#occor.r.group1[occor.p.group1>0.05|abs(occor.r.group1)<0.6] = 0 #remove those r<0.6, p>0.05.
occor.r.group1 <- adj_marine
node_names <- rownames(occor.r.group1)
node_names_clean <- gsub(":.*", "", node_names)
rownames(occor.r.group1) <- node_names_clean
colnames(occor.r.group1) <- node_names_clean
igraph.group1 = graph_from_adjacency_matrix(occor.r.group1,mode="undirected",weighted=TRUE,diag=FALSE)
igraph.group1.weight = E(igraph.group1)$weight
E(igraph.group1)$weight = NA
# set edge color postive correlation pink color, negative blue.
E.color.group1 = igraph.group1.weight
E.color.group1 = ifelse(E.color.group1>0, "pink",ifelse(E.color.group1<0, "blue","grey")) 
E(igraph.group1)$color = as.character(E.color.group1)
#change edge width
E(igraph.group1)$width = abs(igraph.group1.weight)
#change node size if need
#igraph.size.group1 = group1[V(igraph.group1)$name,] 
#igraph.size1.group1 = rowMeans(igraph.size.group1[,2:4])*100
#V(igraph.group1)$size = igraph.size1.group1
# set vertices color, modularity
fc.group1 = cluster_fast_greedy(igraph.group1,weights =NULL)
modularity.group1 = modularity(igraph.group1,membership(fc.group1))
comps.group1 = membership(fc.group1)
colbar.group1 = rainbow(max(comps.group1))
V(igraph.group1)$color = colbar.group1[comps.group1]
#===============group2
#occor.group2 = corr.test(t(group2),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
#occor.r.group2= occor.group2$r 
#occor.p.group2= occor.group2$p 
#occor.r.group2[occor.p.group2>0.05|abs(occor.r.group2)<0.6] = 0
occor.r.group2 <- adj_sediment
node_names <- rownames(occor.r.group2)
node_names_clean <- gsub(":.*", "", node_names)
rownames(occor.r.group2) <- node_names_clean
colnames(occor.r.group2) <- node_names_clean
igraph.group2 = graph_from_adjacency_matrix(occor.r.group2,mode="undirected",weighted=TRUE,diag=FALSE)
igraph.group2.weight = E(igraph.group2)$weight
E(igraph.group2)$weight = NA
E.color.group2 = igraph.group2.weight
E.color.group2 = ifelse(E.color.group2>0, "pink",ifelse(E.color.group2<0, "blue","grey")) 
E(igraph.group2)$color = as.character(E.color.group2)
E(igraph.group2)$width = abs(igraph.group2.weight)
fc.group2 = cluster_fast_greedy(igraph.group2,weights =NULL)
modularity.group2 = modularity(igraph.group2,membership(fc.group2))
comps.group2 = membership(fc.group2)
colbar.group2 = rainbow(max(comps.group2))
V(igraph.group2)$color = colbar.group2[comps.group2]
#===============group3
#occor.group3 = corr.test(t(group3),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
#occor.r.group3= occor.group3$r 
#occor.p.group3= occor.group3$p 
#occor.r.group3[occor.p.group3>0.05|abs(occor.r.group3)<0.6] = 0
occor.r.group3 <- adj_coral
node_names <- rownames(occor.r.group3)
node_names_clean <- gsub(":.*", "", node_names)
rownames(occor.r.group3) <- node_names_clean
colnames(occor.r.group3) <- node_names_clean
igraph.group3 = graph_from_adjacency_matrix(occor.r.group3,mode="undirected",weighted=TRUE,diag=FALSE)
igraph.group3.weight = E(igraph.group3)$weight
E(igraph.group3)$weight = NA
E.color.group3 = igraph.group3.weight
E.color.group3 = ifelse(E.color.group3>0, "pink",ifelse(E.color.group3<0, "blue","grey")) 
E(igraph.group3)$color = as.character(E.color.group3)
E(igraph.group3)$width = abs(igraph.group3.weight)

fc.group3 = cluster_fast_greedy(igraph.group3,weights =NULL)
modularity.group3 = modularity(igraph.group3,membership(fc.group3))
comps.group3 = membership(fc.group3)
colbar.group3 = rainbow(max(comps.group3))
V(igraph.group3)$color = colbar.group3[comps.group3]
#===============group4
#occor.group4 = corr.test(t(group4),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
#occor.r.group4= occor.group4$r 
#occor.p.group4= occor.group4$p 
#occor.r.group4[occor.p.group4>0.05|abs(occor.r.group4)<0.6] = 0
occor.r.group4 <- adj_sponge
node_names <- rownames(occor.r.group4)
node_names_clean <- gsub(":.*", "", node_names)
rownames(occor.r.group4) <- node_names_clean
colnames(occor.r.group4) <- node_names_clean
igraph.group4 = graph_from_adjacency_matrix(occor.r.group4,mode="undirected",weighted=TRUE,diag=FALSE)
igraph.group4.weight = E(igraph.group4)$weight
E(igraph.group4)$weight = NA
E.color.group4 = igraph.group4.weight
E.color.group4 = ifelse(E.color.group4>0, "pink",ifelse(E.color.group4<0, "blue","grey")) 
E(igraph.group4)$color = as.character(E.color.group4)
E(igraph.group4)$width = abs(igraph.group4.weight)

fc.group4 = cluster_fast_greedy(igraph.group4,weights =NULL)
modularity.group4 = modularity(igraph.group4,membership(fc.group4))
comps.group4 = membership(fc.group4)
colbar.group4 = rainbow(max(comps.group4))
V(igraph.group4)$color = colbar.group4[comps.group4]


#################plot the graphs
par(mfrow=c(1,4),mar=c(2,2,2,2))
set.seed(123)
plot(igraph.group1,main="Ambient",vertex.frame.color=NA,vertex.label=NA,edge.width=0.2,
     vertex.size=5, edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
set.seed(123)
plot(igraph.group2,main="Mesophilic Low-solids",vertex.frame.color=NA,vertex.label=NA,edge.width=0.2,
     vertex.size=5, edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))#,layout=layout.graphopt
set.seed(123)
plot(igraph.group3,main="Mesophilic",vertex.frame.color=NA,vertex.label=NA,edge.width=0.2,
     vertex.size=5, edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
set.seed(123)
plot(igraph.group4,main="Mesophilic Co-digest",vertex.frame.color=NA,vertex.label=NA,edge.width=0.2,
     vertex.size=5, edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))#,layout=layout.graphopt

##############random networks
igraph.random.group5=erdos.renyi.game(length(V(igraph.group5)),length(E(igraph.group5)),type="gnm")
igraph.random.group4=erdos.renyi.game(length(V(igraph.group4)),length(E(igraph.group4)),type="gnm")
igraph.random.group3=erdos.renyi.game(length(V(igraph.group3)),length(E(igraph.group3)),type="gnm")
igraph.random.group2=erdos.renyi.game(length(V(igraph.group2)),length(E(igraph.group2)),type="gnm")
igraph.random.group1=erdos.renyi.game(length(V(igraph.group1)),length(E(igraph.group1)),type="gnm")

##############characteristics of network
#group 1 as example
igraph.group1.charac=c( transitivity(igraph.group1), #clustering.coefficient
                        transitivity(igraph.random.group1), #clustering coefficient of random networks
                        modularity(igraph.group1,membership(fc.group1)),#modularity
                        sum(igraph.group1.weight>0)/(sum(igraph.group1.weight>0)+sum(igraph.group1.weight<0)),#positive.ratio
                        mean(igraph::degree(igraph.group1)),#ave.degree
                        average.path.length(igraph.group1), #ave.path.length
                        diameter(igraph.group1, directed = FALSE, unconnected = TRUE, weights = NULL),#diameter
                        length(E(igraph.group1)),#size
                        length(V(igraph.group1)),#order
                        edge_density(igraph.group1,loops=FALSE),#edge.density
                        edge_connectivity(igraph.group1),#edge.connectivity
                        no.clusters(igraph.group1),#number.cluster
                        mean(degree(igraph.group1)/length(V(igraph.group1))),#normalized degree
                        mean(betweenness(igraph.group1,normalized=T)),#normalized betweenness
                        centralization.betweenness(igraph.group1)$centralization) #betweenness.centralization

###############Normalized betweenness and normalized degree of individual genera
degree(igraph.group1)/length(V(igraph.group1))#normalized degree
betweenness(igraph.group1,normalized=T)#normalized betweenness



#  calculate Pi, Zi
library(igraph) 

#igraph.group1
deg_polar <- igraph::degree(igraph.group1)

names(deg_polar) <- sub(":.*", "", names(deg_polar))
seqdeg <- as.numeric(deg_polar)
Nnodes=length(seqdeg)
Z=seqdeg
Z[]=0
P=Z

Membership=membership(fc.group1)
Seq=seq(1:Nnodes)
for(i in 1:Nnodes){
  L=Membership==Membership[i]         
  neighbs=neighbors(igraph.group1,i)               
  Kis=sum(L[neighbs])
  SUM=0
  SUMsq=0	
  SUMP=0
  Miv=Seq[L]
  for(j in 1:sum(L)){
    neighbsj=neighbors(igraph.group1,Miv[j])
    Kjs=sum(L[neighbsj])
    SUM=SUM+Kjs
    SUMsq=SUMsq+Kjs^2
  }
  Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
  if(Kis-SUM/sum(L)==0){Z[i]=0}
  for(k in 1:max(Membership)){
    Lp=Membership==k
    Kisp=sum(Lp[neighbs])
    SUMP=SUMP+(Kisp/seqdeg[i])^2}
  P[i]=1-SUMP
}
attribute_node.group1=cbind(degree=seqdeg,module=Membership,Pi=P,Zi=Z)


#igraph.group2
deg_temperate <- igraph::degree(igraph.group2)

names(deg_temperate) <- sub(":.*", "", names(deg_temperate))
seqdeg <- as.numeric(deg_temperate)
Nnodes=length(seqdeg)
Z=seqdeg
Z[]=0
P=Z

Membership=membership(fc.group2)
Seq=seq(1:Nnodes)
for(i in 1:Nnodes){
  L=Membership==Membership[i]         
  neighbs=neighbors(igraph.group2,i)               
  Kis=sum(L[neighbs])
  SUM=0
  SUMsq=0	
  SUMP=0
  Miv=Seq[L]
  for(j in 1:sum(L)){
    neighbsj=neighbors(igraph.group2,Miv[j])
    Kjs=sum(L[neighbsj])
    SUM=SUM+Kjs
    SUMsq=SUMsq+Kjs^2
  }
  Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
  if(Kis-SUM/sum(L)==0){Z[i]=0}
  for(k in 1:max(Membership)){
    Lp=Membership==k
    Kisp=sum(Lp[neighbs])
    SUMP=SUMP+(Kisp/seqdeg[i])^2}
  P[i]=1-SUMP
}
attribute_node.group2=cbind(degree=seqdeg,module=Membership,Pi=P,Zi=Z)

#igraph.group3
deg_tropical <- igraph::degree(igraph.group3)

names(deg_tropical) <- sub(":.*", "", names(deg_tropical))
seqdeg <- as.numeric(deg_tropical)
Nnodes=length(seqdeg)
Z=seqdeg
Z[]=0
P=Z

Membership=membership(fc.group3)
Seq=seq(1:Nnodes)
for(i in 1:Nnodes){
  L=Membership==Membership[i]         
  neighbs=neighbors(igraph.group3,i)               
  Kis=sum(L[neighbs])
  SUM=0
  SUMsq=0	
  SUMP=0
  Miv=Seq[L]
  for(j in 1:sum(L)){
    neighbsj=neighbors(igraph.group3,Miv[j])
    Kjs=sum(L[neighbsj])
    SUM=SUM+Kjs
    SUMsq=SUMsq+Kjs^2
  }
  Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
  if(Kis-SUM/sum(L)==0){Z[i]=0}
  for(k in 1:max(Membership)){
    Lp=Membership==k
    Kisp=sum(Lp[neighbs])
    SUMP=SUMP+(Kisp/seqdeg[i])^2}
  P[i]=1-SUMP
}
attribute_node.group3=cbind(degree=seqdeg,module=Membership,Pi=P,Zi=Z)


#igraph.group3
deg_sponge <- igraph::degree(igraph.group4)

names(deg_sponge) <- sub(":.*", "", names(deg_sponge))
seqdeg <- as.numeric(deg_sponge)
Nnodes=length(seqdeg)
Z=seqdeg
Z[]=0
P=Z

Membership=membership(fc.group4)
Seq=seq(1:Nnodes)
for(i in 1:Nnodes){
  L=Membership==Membership[i]         
  neighbs=neighbors(igraph.group4,i)               
  Kis=sum(L[neighbs])
  SUM=0
  SUMsq=0	
  SUMP=0
  Miv=Seq[L]
  for(j in 1:sum(L)){
    neighbsj=neighbors(igraph.group4,Miv[j])
    Kjs=sum(L[neighbsj])
    SUM=SUM+Kjs
    SUMsq=SUMsq+Kjs^2
  }
  Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
  if(Kis-SUM/sum(L)==0){Z[i]=0}
  for(k in 1:max(Membership)){
    Lp=Membership==k
    Kisp=sum(Lp[neighbs])
    SUMP=SUMP+(Kisp/seqdeg[i])^2}
  P[i]=1-SUMP
}
attribute_node.group4=cbind(degree=seqdeg,module=Membership,Pi=P,Zi=Z)


#zi pi graph on one plot

png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/Rendered/connectivity_plot_SO.png", width = 5000, height = 3000, res = 600)

par(mfrow=c(1,1),mar=c(4,4,2,14))
plot(attribute_node.group1[,3],attribute_node.group1[,4]
     ,xlim=c(0,1),ylim=c(-4,4),
     xlab="Among-module connectivity P",
     ylab=("Within-module connectivity Z"),
     col=2,
     pch=1,
     cex=0.95,
     main = "")
abline(v=0.62,h=2.5,col=8)
points(attribute_node.group2[,3],attribute_node.group2[,4],
       col=3,pch=6,cex=0.95)
points(attribute_node.group3[,3],attribute_node.group3[,4],
       col=4,pch=2,cex=0.95)
points(attribute_node.group4[,3],attribute_node.group4[,4],
       col=6,pch=3,cex=0.95)
text(0.15,4,"Module hubs")
text(0.8,4,"Network hubs")
text(0.15,-4,"Peripherals")
text(0.8,-4,"Connectors")

legend(1.05, 4,
       legend=c("Marine", "Sediment", "Coral", "Sponge"),
       pch=c(1, 6, 2, 3),  # pch values corresponding to each group
       col=c(2, 3, 4, 6),  # colors corresponding to each group
       xpd=TRUE, 
       bty="n", 
       pt.lwd=2)
dev.off()