library(NST)
# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(vegan)  # For performing PCA
library(dplyr)
library(NetCoMi)
library(microbiome)
library(RColorBrewer)
library(iCAMP)


remove_prefixes <- function(physeq_1) {
  tax_table1 = tax_table(physeq_1)
  taxa_names <- gsub("^d_", "", tax_table1) # Remove 'k_' prefix
  taxa_names <- gsub("^p_", "", taxa_names)
  taxa_names <- gsub("^c_", "", taxa_names)
  taxa_names <- gsub("^o_", "", taxa_names)
  taxa_names <- gsub("^f_", "", taxa_names)
  taxa_names <- gsub("^s_", "", taxa_names) # Remove 's_' prefix
  taxa_names <- gsub("^_", "", taxa_names) # Remove 's_' prefix
  tax_names2 <- tax_table(taxa_names)
  otu_table1 <- otu_table(physeq_1)
  metadata <- sam_data(physeq_1)
  physeq2 <- phyloseq(otu_table1, tax_names2, metadata)
  
  
  return(physeq2)
}
rm(list = setdiff(ls(), c("physeq_combined_Algo")))
gc()  # Run garbage collection to free up memory

combined_data1 <- remove_prefixes(physeq_combined_Algo)
combined_data1 <- physeq_combined_Algo
combined_genus_agglomeration1 <- combined_data1
complete_data1 <- prune_taxa(taxa_sums(combined_genus_agglomeration1) > 0, combined_genus_agglomeration1)

complete_data1 <- core(complete_data1, detection=50, prevalence = 0.01, include.lowest = FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data1 <- prune_samples(non_zero_samples, complete_data1)

combined_N <- complete_data1

otu_mat <- as(otu_table(combined_N), "matrix")
if (taxa_are_rows(combined_N)) {
  otu_mat <- t(otu_mat)
}
tax_table_df <- data.frame(tax_table(combined_N))
meta <- as.data.frame(sample_data(combined_N))
meta$Geographic_Range <- as.factor(meta$Geographic_Range)  # ensure it's a factor

group_df <- data.frame(Geographic_Range = as.factor(meta$Geographic_Range))

rownames(group_df) <- rownames(meta)  # Important!

nst_result <- tNST(
  comm = otu_mat,
  group = group_df,
  dist.method = "bray",
  rand = 200,
  meta.com=NULL, abundance.weighted=TRUE,
  output.rand=TRUE, LB=FALSE, null.model="PF",
  between.group=FALSE, SES=TRUE, RC=TRUE
)

saveRDS(nst_result, "../nst_result.rds")
nst_result = readRDS("../combined_N_NST.rds")
tnstbt=nst.boot(nst.result=nst_result, group=group_df, rand=500, trace=TRUE,
                two.tail=FALSE, out.detail=TRUE, between.group=FALSE)

saveRDS(tntbt, "../tnstsbt.rds")

tnstbt <- readRDS("../tnstsbt.rds")

nst_data <- tnstbt$summary[tnstbt$summary$Index == "NST", ]
ggplot(nst_data, aes(x = Group, ymin = LowerWhisker, lower = LowerHinge, 
                     middle = Median, upper = HigherHinge, ymax = HigherWhisker)) +
  geom_boxplot(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "NST Boxplot for Polar, Temperate, and Tropical",
       y = "Value", x = "Climate Zone") +
  theme_minimal()

