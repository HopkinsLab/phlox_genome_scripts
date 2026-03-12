#!/usr/bin/env Rscript
rm(list=ls())

# library(R.utils)
library(data.table)
library(stringr)

argv <- commandArgs(trailingOnly=TRUE)
# argv <- c("results/multi/cusp_drum_roem_pilo/res_stacks/min_samples_per_pop~0.5/min_samples_overall~0.5/mind_0.75/populations.all.lFilt.iFilt.pop", "1")

pop_fn <- argv[1]
out_fn <- argv[2]

admixed_samps <- c("640-2-1", "651-8-1", "647-5-1", "666-6-1", "640-2-1", "725_12", "678-3-1")
north_pilo_samps <- c("P_pilosa_pilosa_OPGC-4183", "P_pilosa_pilosa_PIL-03-01", "P_pilosa_pilosa_PIL_06_1983", "P_pilosa_pilosa_PIL_09_1929", "P_pilosa_pilosa_PIL_10_1955", "P_pilosa_pilosa_PIL_11_1924")

rm_samps <- c(admixed_samps, north_pilo_samps)

popd <- fread(pop_fn, header=FALSE, sep="\t", col.names=c("ID", "SPECIES"))
popd <- popd[!(ID %in% rm_samps)] # Filter admixed and northern-pilosa samples
setkey(popd, SPECIES)

outgrp <- popd["pilo", ID] # Get outgroup individiuals
if("amoe" %in% popd$SPECIES){
    outgrp <- c(outgrp, popd["amoe", ID])
}

trio_sp <- c("roem", "drum", "cusp")

# Load data
dat <- fread("file:///dev/stdin", sep="\t", header=TRUE)

for(samp in rm_samps){
    if(samp %in% names(dat)){
        dat[, eval(samp) := NULL]
    }
}

# Filter to keep only sites where all pilosa samples are ref
mask <- dat[, outgrp, with=FALSE][, all(unlist(.SD) %in% c("0", "./.")), by=1:nrow(dat)]$V1
dat <- dat[mask]; gc()
tmp <- dat[, outgrp, with=FALSE][, sum(unlist(.SD) == "0"), by=1:nrow(dat)]$V1
# dat[, Outgroup_N_GT := tmp]
# dat <- dat[Outgroup_N_GT > 0]; gc()
dat <- dat[tmp > 0]; gc()


# For each drummondii sample, choose a random roem and cusp to compare against (3X)
trio_m <- list()
for(rsamp in popd["roem",ID]){
    for(dsamp in popd["drum",ID]){
        trio_m[[length(trio_m)+1]] <- data.table(roem=rsamp, drum=dsamp, cusp=popd["cusp",ID])
    }
}
trio_m <- rbindlist(trio_m)

# Keep only samples we're going to be analyzing
samps <- unique(as.vector(as.matrix(trio_m)))
gts <- as.matrix(dat[, samps, with=FALSE])
# dat <- dat[, .(CHROM, POS, REF, ALT, Outgroup_N_GT)]; gc()
rm(dat); gc()

# Function for computing patterns
CountPatterns <- function(samps){
    m <- gts[,samps]
    
    for(i in 1:3){
        m <- m[m[,i] != "./.",]
    }
    
    cur_mask <- m %in% c("0/1", "1/0")
    nhet <- sum(cur_mask)
    m[cur_mask] <- sample(as.character(0:1), size=nhet, replace=TRUE)
    
    counts <- table(factor(as.numeric(m[,3]) + 2*as.numeric(m[,2]) + 4*as.numeric(m[,1]), levels=0:7))
    rm(m); gc()
    # return(counts)
    return(as.list(counts))
}
patt_counts <- trio_m[, CountPatterns(unlist(.SD)), by=1:nrow(trio_m)]
rm(gts); gc()
patt_counts[, nrow := NULL]
# patt_counts <- apply(trio_m, 1, CountPatterns)

# Write out counts
DecToBinPatt <- function(x){
    return( str_c(c(rev(as.numeric(intToBits(x))[1:3]), "0"), collapse="") )
}

names(patt_counts) <- sapply(0L:7L, DecToBinPatt)
for(nm in names(trio_m)){
    patt_counts[, eval(nm) := trio_m[[nm]]]
}

fwrite(patt_counts, file=out_fn, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
