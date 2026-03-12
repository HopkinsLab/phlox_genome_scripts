#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)

argv <- commandArgs(trailingOnly=TRUE)

frq_d <- fread(argv[1], header=TRUE)

frq_d[, ID := str_c(CHR, SNP, sep="_")]
frq_d[, AC_STR := str_c(MAC, NCHROBS-MAC, sep=",")]

out_d <- dcast(frq_d, ID ~ CLST, value.var="AC_STR")
out_d[, ID := NULL]

fwrite(out_d, file="/dev/stdout", sep=" ", col.names=TRUE, row.names=FALSE, quote=FALSE)
