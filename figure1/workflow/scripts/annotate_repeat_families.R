#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)

argv <- commandArgs(trailingOnly = TRUE)

bed_fn <- argv[1] #"results/annotations/repeats.bed"
rpt_bed <- fread(bed_fn, header=FALSE)
cnames <- c("CHROM", "START", "END", "NAME", "SCORE", "STRAND", "V7", "V8", "V9", "INFO")
names(rpt_bed) <- cnames
rpt_bed[, ID := str_remove(str_remove(str_extract(INFO, '"Motif:.*"'), '"Motif:'), '"')]

fai_fn <- argv[2] #"/n/holylfs05/LABS/hopkins_lab/Lab/Phlox_resources/WG_assemblies/annotations/drummondii/phlox_flye.v1.0_combined_repeats.fa.fai"
rpt_fai <- fread(cmd=str_interp("cut -f 1 ${fai_fn}"), sep="#", header=FALSE)
names(rpt_fai) <- c("ID", "FAMILY")
rpt_fai <- rpt_fai[ID != "pdrum_sat"]

rpt_fai[str_detect(FAMILY, "/"), c("SUPERFAM", "SUBFAM") := tstrsplit(FAMILY, "/")]

rpt_bed <- merge(rpt_bed, rpt_fai, by="ID", all.x=TRUE)

# rpt_bed[!str_detect(ID, "family"), unique(ID)]
rpt_bed[str_detect(ID, "^\\([ATCG]+\\)n"), SUPERFAM := "Simple_repeat"]
rpt_bed[str_detect(ID, "-rich"), SUPERFAM := "Low_complexity"]
rpt_bed[ID == "polypurine", SUPERFAM := "Low_complexity"]
rpt_bed[ID == "pdrum_sat", SUPERFAM := "Satellite"]

x <- c("Unknown", "rRNA", "LTR", "Satellite", "DNA", "snRNA", "Simple_repeat", "tRNA")
rpt_bed[FAMILY %in% x, SUPERFAM := FAMILY]

out_fn <- str_replace(bed_fn, "\\.bed$", ".families.bed")
setcolorder(rpt_bed, c(cnames, "FAMILY", "SUPERFAM", "SUBFAM", "ID"))
fwrite(
    rpt_bed, 
    file=out_fn, 
    col.names=FALSE, 
    row.names=FALSE, 
    sep="\t", 
    quote=FALSE, 
    na=".")
