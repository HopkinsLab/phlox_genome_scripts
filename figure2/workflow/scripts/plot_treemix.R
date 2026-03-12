#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(R.utils)
library(stringr)

argv <- commandArgs(trailingOnly=TRUE)

src_dir <- argv[1]
prefix <- argv[2]

src_dir <- "/n/holylfs05/LABS/hopkins_lab/Lab/software/treemix-1.13/src"
prefix <- "results/multi/cusp_drum_roem_pilo/res_stacks/min_samples_per_pop~0.1/min_samples_overall~0.1/mind_0.9/treemix/populations.snps.lFilt.iFilt.input"
tmpdir <- "results/tmp"

source(str_interp("${src_dir}/plotting_funcs.R"))

# Plot tree
out_fn <- str_interp("${prefix}.tree.pdf")
pdf(file=out_fn, width=8, height=4)
tree_dat <- plot_tree(prefix)
dev.off()
print(out_fn)

# # Plot residuals
# pops <- unlist(fread(str_interp("${prefix}.gz"),  nrows=1, header=FALSE)[1,])
# poporder_fn <- str_interp("${tmpdir}/${sample(10000,1)}.poporder")
# fwrite(data.table(V1=pops), file=poporder_fn, col.names=FALSE, row.names=FALSE, quote=FALSE)

# out_fn <- str_interp("${prefix}.resid.pdf")
# pdf(file=out_fn, width=8, height=8)
# resid_dat <- plot_resid(prefix, poporder_fn)
# dev.off()
