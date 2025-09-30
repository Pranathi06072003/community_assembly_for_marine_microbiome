# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(vegan)  # For performing PCA
library(dplyr)
library(ape)
library(NetCoMi)
library(microbiome)
remove_prefixes <- function(physeq_1) {
  tax_table1 = tax_table(physeq_1)
  taxa_names <- gsub("^d_|^p_|^c_|^o_|^f_|^s_|^_", "", tax_table1)
  tax_names2 <- tax_table(taxa_names)
  otu_table1 <- otu_table(physeq_1)
  metadata <- sample_data(physeq_1)
  tree <- phy_tree(physeq_1)  # preserve the tree
  
  physeq2 <- phyloseq(otu_table1, tax_names2, metadata, tree)  # include tree
  
  return(physeq2)
}

ps_temperate <- subset_samples(physeq_entire, Geographic_Range == "Temperate")
ps_tropical <- subset_samples(physeq_entire, Geographic_Range == "Tropical")
ps_polar <- subset_samples(physeq_entire, Geographic_Range == "Polar")

rm(list = setdiff(ls(), c("ps_temperate", "ps_tropical", "ps_polar")))
gc()  # Run garbage collection to free up memory

ps_polar <- core(ps_polar, detection=50, prevalence = 0.01, include.lowest = FALSE)
ps_tropical <- core(ps_tropical, detection=50, prevalence = 0.1, include.lowest = FALSE)
ps_temperate <- core(ps_temperate, detection=50, prevalence = 0.1, include.lowest = FALSE)

non_zero_samples <- sample_sums(ps_temperate) > 0
complete_data1 <- prune_samples(non_zero_samples, ps_temperate)
complete_data1 <- prune_taxa(taxa_sums(complete_data1) > 0, complete_data1)

combined_N <- transform_sample_counts(complete_data1, function(x) x / sum(x))

combined_data <- complete_data1

phy_tree <- phy_tree(complete_data1)
otu_table <- otu_table(complete_data1)
tax_table <- tax_table(complete_data1)
env_data <- sample_data(complete_data1)
env_df <- data.frame(env_data)

View(env_df)
treat <- env_df[, "sample_origin_code", drop = FALSE]

if (taxa_are_rows(complete_data1)) {
  otu_table <- t(otu_table)
}
write.csv(otu_table, file = "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/otu_table_icamp.csv", quote = FALSE)
write.csv(tax_table, file = "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/tax_table_icamp.csv", quote = FALSE)
write.tree(phy_tree, file = "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/phylogenetic_tree_icamp.nwk")
write.csv(env_df, file = "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/env_table_icamp.csv", quote = FALSE)

save.wd="C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/wd_icamp8"
if(!dir.exists(save.wd)){dir.create(save.wd)}

rand.time=100  # randomization time, 1000 is usually enough. For example test, you may set as 100 or less to save time.
nworker=1 # nworker is thread number for parallel computing, which depends on the CPU core number of your computer.
memory.G=16 # to set the memory size as you need (but should be less than the available space in your hard disk), so that calculation of large tree will not be limited by physical memory. unit is Gb.

# 3 # load R packages and data
library(iCAMP)
library(ape)

comm <- read.table("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/otu_table_icamp.csv",
                   header = TRUE,
                   sep = ",",
                   row.names = 1,
                   as.is = TRUE,
                   stringsAsFactors = FALSE,
                   comment.char = "",
                   check.names = FALSE)
tree=read.tree(file = "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/phylogenetic_tree_icamp.nwk")

clas <- read.table("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/tax_table_icamp.csv", header = TRUE, sep = ",", row.names = 1,
                   as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                   check.names = FALSE)

env <- read.table("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/env_table_icamp.csv", header = TRUE, sep = ",", row.names = 1,
                  stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE)

sampid.check=match.name(rn.list=list(comm=comm,treat=treat,env=env))
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE] # if some unmatched samples were removed, some OTUs may become ghosts, then you may use this line to remove them if necessary.
env=sampid.check$env # skip this if you do not have env.file

spid.check=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas),tree.list=list(tree=tree))

comm=spid.check$comm
clas=spid.check$clas
tree=spid.check$tree

# 6 # calculate pairwise phylogenetic distance matrix.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large phylogenetic distance matrix. 
setwd(save.wd)
if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
  # output files:
  # path.rda: a R object to list all the nodes and  edge lengthes from root to every tip. saved in R data format. an intermediate output when claculating phylogenetic distance matrix.
  # pd.bin: BIN file (backingfile) generated by function big.matrix in R package bigmemory. This is the big matrix storing pairwise phylogenetic distance values. By using this bigmemory format file, we will not need memory but hard disk when calling big matrix for calculation.
  # pd.desc: the DESC file (descriptorfile) to hold the backingfile (pd.bin) description.
  # pd.taxon.name.csv: comma delimited csv file storing the IDs of tree tips (OTUs), serving as the row/column names of the big phylogenetic distance matrix.
}else{
  # if you already calculated the phylogenetic distance matrix in a previous run
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}

# you may skip step 7-8, if the "alternative way" based on stochasticity is applicable, as mentioned in the method part of iCAMP paper (Ning et al 2020 Nature Communications).
# 7 # assess niche preference difference between species
# env is required for this step.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large niche difference matrix. 
setwd(save.wd)

env_numeric <- env
env_numeric <- data.frame(env[, c("latitude_deg", "longitude_deg")], 
                          stringsAsFactors = FALSE)

env_numeric$latitude_deg <- as.numeric(env_numeric$latitude_deg)
env_numeric$longitude_deg <- as.numeric(env_numeric$longitude_deg)
env_use <- na.omit(env_numeric)

# Convert to matrix for iCAMP
env_use <- as.matrix(env_use)

# Check final result
str(env_use)


niche.dif <- iCAMP::dniche(env = env_use, comm = comm, method = "niche.value",
                           nworker = 3, out.dist = FALSE, bigmemo = TRUE,
                           nd.wd = save.wd)


# 8.1 # try phylogenetic binning using current setttings.
ds = 0.2 # setting can be changed to explore the best choice
bin.size.limit = 48 # setting can be changed to explore the best choice. # here set as 5 just for the small example dataset. For real data, usually try 12 to 48.

# The tree for taxa.binphy.big must be a rooted tree.
if(!ape::is.rooted(tree))
{
  tree.rt=iCAMP::midpoint.root.big(tree = tree, pd.desc = pd.big$pd.file,
                                   pd.spname = pd.big$tip.label,pd.wd = pd.big$pd.wd,
                                   nworker = nworker)
  tree=tree.rt$tree
}
phylobin=taxa.binphy.big(tree = tree, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = 3)

sp.bin=phylobin$sp.bin[,3,drop=FALSE]
sp.ra=colMeans(comm/rowSums(comm))
abcut=3 # you may remove some species, if they are too rare to perform reliable correlation test.
commc=comm[,colSums(comm)>=abcut,drop=FALSE]
dim(commc)
spname.use=colnames(commc)

env <- env_use
binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "pearson",r.cut = 0.1, p.cut = 0.05, min.spn = 5)

sig_bins <- binps$detail[
  binps$detail[,"latitude_deg.pearson.p"] < 0.05 |
    binps$detail[,"longitude_deg.pearson.p"] < 0.05,
]
sig_abundant_bins <- sig_bins[sig_bins[,"abu"] > 0.01, ]


# 9 # iCAMP analysis
# 9.1 # without omitting small bins.
# commonly use # set sig.index as Confidence instead of SES.RC (betaNRI/NTI + RCbray)
bin.size.limit = 24 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
sig.index="Confidence" # see other options in help document of icamp.big.
prefix="Test"  # prefix of the output file names. usually use a project ID.
icres_temperate=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                                 pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                                 prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                                 phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                                 phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                                 nworker = 3, memory.G = 15, rtree.save = FALSE, detail.save = TRUE, 
                                 qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                                 correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                                 ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
#
temp_icamp_result <- icres_temperate$detail$processes$CbMPDiCBraya

treat<- as.vector(treat)
levels(factor(treat))
icbin_temperate=icamp.bins(icamp.detail = icres_temperate$detail,treat=NULL
                           clas = clas,boot = TRUE, rand.time = 200, between.group = FALSE)

process_values <- icbin_temperate$Pt[, c("HeS", "HoS", "DL", "HD", "DR")]
process_vector <- as.numeric(process_values)
process_df <- data.frame(Process = c("HeS", "HoS", "DL", "HD", "DR"), Value = process_vector)

process_df$Process <- factor(process_df$Process, levels = c("DL", "DR", "HD", "HoS", "HeS"))

# Filter out rows with 0 value
process_df <- subset(process_df, Value > 0)

# Reorder the factor levels based on what's left
process_df$Process <- factor(process_df$Process, levels = c("DL", "DR", "HD", "HoS", "HeS"))

# Sort the dataframe
process_df <- process_df[order(process_df$Process), ]

# Compute fractions and label positions
process_df$Fraction <- process_df$Value / sum(process_df$Value)
process_df$ymax <- cumsum(process_df$Fraction)
process_df$ymin <- c(0, head(process_df$ymax, n = -1))
process_df$LabelPosition <- (process_df$ymin + process_df$ymax) / 2
process_df$PercentLabel <- paste0(round(process_df$Fraction * 100, 1), "%")

# Load ggplot2
library(ggplot2)

# Plot
ggplot(process_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Process)) +
  geom_rect(color = "white") +
  coord_polar(theta = "y") +
  xlim(c(2, 4.5)) +
  theme_void() +
  geom_text(aes(x = 4.3, y = LabelPosition, label = PercentLabel), size = 4) +
  scale_fill_brewer(palette = "Set2", drop = TRUE) +  # drop unused levels from legend
  ggtitle("") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

# Save the plot
ggsave("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/QZA_Results/temperate_icamp.png",
       width = 10, height = 6, dpi = 800)
