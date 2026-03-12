#!/usr/bin/env Rscript
# make a "Super-population" popmap file based on lat-long coords

library(data.table)
library(stringr)

argv <- commandArgs(trailingOnly=TRUE)

samps <- fread(argv[1], header=FALSE, sep=" ")$V1 #pop file
met <- fread(argv[2], header=TRUE, sep="\t")

setkey(met, Sample_ID)
met <- met[(samps)]

pops <- met[, .(lat=lat[1], long=long[1], species=species[1]), by=Population]
setkey(pops, Population)

# Get super populations
cutoff <- 1e-4
sup_pops <- list()
obs_pops <- c()
for(p in pops[, Population]){
    if(p %in% obs_pops){
        next
    }
    ll <- pops[(p), c(lat, long)]
    sp <- pops[(p), species]
    same_p <- pops[(abs(lat - ll[1]) < cutoff) & (abs(long - ll[2]) < cutoff) & (species == sp), unique(Population)]
    
    sup_pops[[str_c(same_p, collapse="_")]] <- same_p
    obs_pops <- c(obs_pops, same_p)
}

# Write out super populations
out <- met[, .(Sample_ID, Population)]
for(sp in names(sup_pops)){
    out[Population %in% sup_pops[[sp]], SuperPop := sp]
}

fwrite(out[, .(Sample_ID, SuperPop)], file="/dev/stdout", sep="\t", 
       quote=FALSE, row.names=FALSE, col.names=FALSE)
