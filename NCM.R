# --- Required libraries ---
library(phyloseq)
library(vegan)
library(ggplot2)
library(MicEco)
library(Hmisc)
library(minpack.lm)
library(stats4)
library(grid)

# --- Utility functions ---
remove_prefixes <- function(physeq_1) {
  tax_table1 <- tax_table(physeq_1)
  cleaned <- apply(tax_table1, 2, function(x) gsub("^[a-z]__*", "", x))
  tax_table(physeq_1) <- tax_table(cleaned)
  return(physeq_1)
}

fit_neutral_model <- function(ps, out_file) {
  otu_mat <- as.matrix(otu_table(ps))
  spp <- t(otu_mat)
  
  N <- mean(rowSums(spp))
  p.m <- colMeans(spp)
  p.m <- p.m[p.m != 0]
  p <- p.m/N
  
  spp.bi <- 1 * (spp > 0)
  freq <- colMeans(spp.bi)
  freq <- freq[freq != 0]
  
  C <- merge(p, freq, by = 0)
  C <- C[order(C[,2]), ]
  C <- as.data.frame(C)
  C.0 <- C[rowSums(C == 0) == 0, ]
  
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  d <- 1/N
  
  # Fit neutral model
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 - p), lower.tail = FALSE),
                 start = list(m = 0.00001))
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 - p), lower.tail = FALSE)
  pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05,
                     method = "wilson", return.df = TRUE)
  
  Rsqr <- 1 - sum((freq - freq.pred)^2) / sum((freq - mean(freq))^2)
  
  bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[,2:3])
  inter.col <- rep('black', nrow(bacnlsALL))
  inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- '#A52A2A'
  inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- '#29A6A6'
  
  # --- Plotting ---
  png(out_file, width = 3500, height = 3500, res = 800)
  grid.newpage()
  pushViewport(viewport(h=0.6,w=0.6))
  pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
  grid.rect()
  grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch=20, gp=gpar(col=inter.col,cex=0.7))
  grid.yaxis(); grid.xaxis()
  grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp=gpar(col='blue',lwd=2), default='native')
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp=gpar(col='blue',lwd=2,lty=2), default='native') 
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp=gpar(col='blue',lwd=2,lty=2), default='native')  
  grid.text(y=unit(0,'npc')-unit(2.5,'lines'), label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
  grid.text(x=unit(0,'npc')-unit(3,'lines'), label='Frequency of Occurrence', gp=gpar(fontface=2), rot=90) 
  grid.text(paste("R² =", round(Rsqr,3), "\n", "Nm =", round(coef(m.fit)*N)),
            x=unit(0.8,"npc"), y=unit(0.8,"npc"))
  dev.off()
}

# --- Step 1: Preprocess main dataset ---
physeq_clean <- remove_prefixes(physeq_combined_Algo)
physeq_clean <- subset_samples(physeq_clean, sample_origin == "marine metagenome")

# --- Step 2: Run models ---
# Whole dataset
ps_whole <- core(prune_samples(sample_sums(physeq_clean) > 0, physeq_clean),
                 detection = 50, prevalence = 0.001)
fit_neutral_model(ps_whole, "neutral_model_fit_whole.png")

# Tropical
ps_trop <- subset_samples(physeq_clean, Geographic_Range == "Tropical")
ps_trop <- core(prune_samples(sample_sums(ps_trop) > 0, ps_trop),
                detection = 50, prevalence = 0.001)
fit_neutral_model(ps_trop, "neutral_model_fit_tropical.png")

# Temperate
ps_temp <- subset_samples(physeq_clean, Geographic_Range == "Temperate")
ps_temp <- core(prune_samples(sample_sums(ps_temp) > 0, ps_temp),
                detection = 50, prevalence = 0.001)
fit_neutral_model(ps_temp, "neutral_model_fit_temperate.png")

# Polar
ps_polar <- subset_samples(physeq_clean, Geographic_Range == "Polar")
ps_polar <- core(prune_samples(sample_sums(ps_polar) > 0, ps_polar),
                 detection = 50, prevalence = 0.001)
fit_neutral_model(ps_polar, "neutral_model_fit_polar.png")
