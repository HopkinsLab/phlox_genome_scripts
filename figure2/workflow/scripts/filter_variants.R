#!/usr/bin/env Rscript

library(data.table)
library(stringr)

argv <- commandArgs(trailingOnly=TRUE)

min_gq <- argv[1]
min_dp <- argv[2]
max_dp <- argv[3]

# Main snp data (vcf format)
vcf_colnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
dat <- fread("file:///dev/stdin", header=FALSE, sep="\t")
samp_cols <- str_c("S", 1:(ncol(dat) - length(vcf_colnames)))
names(dat) <- c(vcf_colnames, samp_cols)


#############################################
# Filtering

# Parse format string
fmt_names <- strsplit(dat[FORMAT != "GT", FORMAT[1]], ":")[[1]]

# For each sample
for(s in samp_cols){
    mask <- !(dat[[s]] %in% c("0/0", "./.")) # Focus only on non-ref non-missing gts

    cur_stats <- tstrsplit(dat[mask,][[s]], ":")
    names(cur_stats) <- fmt_names

    # Filter on quality stats
    mask <- 
    cur_stats[["GQ"]] <- as.numeric(cur_stats[["GQ"]])
}


dat[, (fmt_names) := tstrsplit(SAMPLE, ":")]
dat[, GQ  := as.integer(GQ)]
dat[, DP  := as.integer(DP)]
dat[, AF  := as.numeric(AF)]
dat[, MAC := as.integer(DP * AF)]

FilterOnQualityStats <- function(in_dat, gt){
    keep <- rep(TRUE, nrow(in_dat))
    for(j in setdiff(names(cutoffs), "GT")){ # For each quality statistic, j
        coff <- cutoffs[gt][[j]]
        if(!is.na(coff)){
            if(j == "AD"){
                mask <- apply(in_dat[, .(MAC, DP)], 1, function(x) binom.test(x[1], x[2], p=0.05)$p.value > coff)
            } else {
                mask <- (in_dat[[j]] > coff)
            }
            keep <- keep & mask
        }
    }
    return(keep)
}

dat[, KEEP := TRUE]
for(gt in cutoffs[["GT"]]){ # For each genotype, gt
    dat[, KEEP := FilterOnQualityStats(.SD, GT[1]), by=GT]
}


# Write out
fwrite(dat[(KEEP)][, vcf_colnames, with=FALSE], 
       file="/dev/stdout", quote=FALSE, sep="\t", 
       row.names=FALSE, col.names=FALSE)
