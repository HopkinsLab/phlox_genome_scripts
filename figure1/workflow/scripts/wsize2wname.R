#!/usr/bin/env Rscript

library(stringr)

argv <- commandArgs(trailingOnly=TRUE)
wname <- argv[1]

N <- str_extract(wname, "^\\d+")
U <- substr(wname, start=str_length(N)+1, stop=str_length(wname))
N <- as.numeric(N)

if(U %in% c("b", "bp")){
    pwr <- 0
} else if(U %in% c("kb", "kbp")){
    pwr <- 3
} else if(U %in% c("Mb", "Mbp")){
    pwr <- 6
} else if(U %in% c("Gb", "Gbp")){
    pwr <- 9
}

out_val <- N * (10**pwr)
cat(str_interp("$[d]{out_val}\n"))
