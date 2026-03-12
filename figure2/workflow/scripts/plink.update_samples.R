#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)

argv <- commandArgs(trailingOnly=TRUE)

met <- fread(argv[1], header=TRUE, sep="\t")
samps <- fread(argv[2], header=FALSE, sep="\t")$V1

setkey(met, Sample_ID)

out_d <- met[(samps), .(Sample_ID, Sample_ID, Population, Sample_ID)]

fwrite(out_d, file="/dev/stdout", sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
