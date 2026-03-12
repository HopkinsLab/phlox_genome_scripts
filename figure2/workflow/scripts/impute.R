#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)
library(ggpubr)
library(missForest)

# Read in args 
argv <- commandArgs(trailingOnly=TRUE)

if(length(argv) == 0){
    argv <- c("results/tmp/example", "results/tmp/example", "10", "8", "8888")
}

in_prefix <- argv[1]
out_prefix <- argv[2]
maxiter <- as.integer(argv[3])
ncores <- as.integer(argv[4])
use_seed <- as.numeric(argv[5])

# Set up parallel compute if need be
if(ncores > 1 ){    
    library(doParallel)
    library(doRNG)
    
    registerDoParallel(cores=ncores)
    if(getDoParWorkers() != ncores){
       stop(str_interp("Tried to register ${ncores} cores, but failed somehow"))
    }

    registerDoRNG(seed = use_seed)

    parallel_arg <- "forests"
} else {
    set.seed(use_seed)
    parallel_arg <- "no"
}

# Read in data
cat(date(), ": Reading in data\n")
in_fn <- str_interp("${in_prefix}.traw")
cnames <- unname(unlist(fread(in_fn, header=FALSE, nrows=1)))
cnames[3] <- "CM"
cnames[-(1:6)] <- unlist(lapply(strsplit(cnames[-(1:6)], "_"), function(x) str_c(x[-1], collapse="_")))
dat <- fread(in_fn, header=FALSE, na.strings="NA", skip=1, col.names=cnames)

# Make sure all columns are categorical
for(nm in names(dat)[-(1:6)]){
    levs <- unique(dat[[nm]])
    dat[, eval(nm) := factor(get(nm), levels=levs[!is.na(levs)])]
}


# Compute priors based on genotype freq
computePrior <- function(x){
    x <- x[!is.na(x)]
    denom <- length(x)

    out_vec <- c()
    for(i in levels(x)){
        out_vec[i] <- sum(x == i) / denom
    }

    return(out_vec)
}

priors <- lapply(dat[,-(1:6)], computePrior)

ntree <- 100
mtry <- floor(sqrt(ncol(dat)))

cat(date(), ": Run missForest\n")
imp <- missForest(dat[,-(1:6)], classwt=priors, maxiter=maxiter, ntree=ntree, mtry=mtry, 
                    parallelize=parallel_arg, verbose=TRUE)

# Write R objects
cat(date(), ": Write R objects\n")

out_fn <- str_interp("${out_prefix}.Rdata")
cat(str_interp("Saving R objects to ${out_fn}\n"))
save(list=ls(), file=out_fn) 

# Confirm that tfam has same sample order as traw
famsamps <- fread(str_interp("${in_prefix}.tfam"), header=FALSE)$V2
if(!all(names(dat)[-(1:6)] == famsamps)){
    stop(str_interp("Samples names in ${in_prefix}.tfam\ndo not match those in ${in_prefix}.traw\n"))
}

# Copy old tfam
system(str_interp("cp ${in_prefix}.tfam ${out_prefix}.tfam"))

# Write new tped
out_dat <- dat[, .(CHR, SNP, CM, POS)]
for(samp in famsamps){
    newcol1 <- str_c(samp, "_1")
    newcol2 <- str_c(samp, "_2")
    
    gts <- imp$ximp[[samp]]
    
    # Init alleles
    a1 <- rep(0L, length(gts))
    a2 <- rep(0L, length(gts))
    
    # Homozygous major (allele code: 2)
    a1[gts == 0] <- 2L
    a2[gts == 0] <- 2L

    # Homozygous minor (allele code: 1)
    a1[gts == 2] <- 1L
    a2[gts == 2] <- 1L

    # Hets (1/2)
    a1[gts == 1] <- 1L
    a2[gts == 1] <- 2L

    # Update
    out_dat[, eval(newcol1) := a1]
    out_dat[, eval(newcol2) := a2]
}

out_fn <- str_interp("${out_prefix}.tped")
fwrite(out_dat, file=out_fn, sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)

# Output error
out_fn <- str_interp("${out_prefix}.OOBerror.txt")
cat(str_interp("$[0.8f]{imp$OOBerror}\n"), file=out_fn)

