#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(viridis)


argv <- commandArgs(trailingOnly = TRUE)
in_prefix <- argv[1]
met_fn <- argv[2]
out_prefix <- argv[3]

max_pcs <- 3 # How many PCs to compare up to

# in_prefix <- "merged.drum_cusp_roem.minGQ_30.pca"
# met_fn <- "metadata.tsv"

# Read in metadata
met_d <- fread(met_fn, header=TRUE, sep="\t", na.strings=c("NA", "#N/A"))
setnames(met_d, names(met_d), str_replace(names(met_d), " ", "_"))

# Read in eigenvectors
evec_d <- fread(str_interp("${in_prefix}.eigenvec"), header=FALSE, sep=" ")
names(evec_d)[1:2] <- c("FID", "IID")
names(evec_d)[3:ncol(evec_d)] <- str_c("PC", 1:(ncol(evec_d) - 2))

# From eigenvalues compute proportion variance explained
eval_d <- fread(str_interp("${in_prefix}.eigenval"), header=FALSE)$V1
prop_var <- eval_d / sum(eval_d)

# Make plot dat
plt_d <- merge(evec_d, met_d, by.x="IID", by.y="Sample_ID")

# Which colors to try
color_groups <- c("species", "region", "lat", "long")

# PCs to plot
pc_names <- str_c("PC", 1:max_pcs)
pc_pairs <- combn(pc_names, 2) # Plot all pairwise comparisons

for(i in 1:ncol(pc_pairs)){
  pcA <- pc_pairs[1,i]
  pcB <- pc_pairs[2,i]

  nA <- as.integer(str_remove(pcA, "PC"))
  nB <- as.integer(str_remove(pcB, "PC"))

  for(color_by in color_groups){
    if(typeof(plt_d[[color_by]]) != "double"){
      plt_d[is.na(get(color_by)), eval(color_by) := "unknown"]
    }
    
    p <- ggplot(plt_d, aes_string(x=pcA, y=pcB)) + 
      geom_point(aes_string(color=color_by) , alpha=0.4) + 
      theme_bw() + theme(legend.position="top") +
      xlab(str_interp("${pcA} ($[0.1f]{prop_var[nA]*100}%)")) + 
      ylab(str_interp("${pcB} ($[0.1f]{prop_var[nB]*100}%)"))
    
    if(typeof(plt_d[[color_by]]) == "double"){
      p <- p + scale_color_viridis_c()
    } else {
      p <- set_palette(p, "d3")
    }
    
    out_fn <- str_interp("${out_prefix}.${pcA}_${pcB}.cBy_${color_by}.pdf")
    ggsave(p, file=out_fn, width=5, height=5)
  }
}
