#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)

argv <- commandArgs(trailingOnly=TRUE)
aggr_method <- tolower(argv[1])

out_dir <- argv[2]

if(aggr_method == "pressubsamp" || aggr_method == "pressub" || aggr_method == "pressamp"){
    if(length(argv) > 2){
        rep_i <- argv[3]
        use_seed <- as.integer(argv[4])
        set.seed(use_seed)
    } else {
        stop("Subsampling requires a replicate id (arg 3) and a random seed (arg 4)\n")
    }
}

# Load data
repeat_summary <- fread('tables/REPEAT_SUMMARY.tab', sep='\t')
unique_repIDs <- unique(repeat_summary[, repID])

sv_summary <- fread('tables/phlox_variation_summary.FINAL.tab.gz', sep='\t')

# we don't need the svs with no repeats associated
sv_summary <- sv_summary[!is.na(repID),]

# we also don't need svs with missing data for sv frequency
sv_summary <- sv_summary[!(AF_pdrum == '.'),][!(AF_proe == '.'),][!(AF_pcusp == '.'),]

CountUnique <- function(x, lev){
    # Count unique occurrences in "x", assuming levels "lev"
    y <- as.vector(table(factor(unique(x[!is.na(x)]), levels=lev)))
    return(data.table(LEVELS=lev, COUNT=y))
}

CountRepIDs <- function(in_dat, lev){
    # Count repeat id occurrences assuming levels "lev"
    vals <- in_dat[, CountUnique(str_split(repID, ",")[[1]], lev), by=ID][, .(COUNT=sum(COUNT)), by=LEVELS]
    out_vec <- vals$COUNT
    names(out_vec) <- vals$LEVELS
    return( out_vec )
}

if(aggr_method == "uniq" || aggr_method == "unique" || aggr_method == "onesp" || aggr_method == "1sp"){
    # Unique (1sp counts)
    out_d <- list()
    out_d[["drumUniq"]] <- CountRepIDs(sv_summary[AF_pdrum > 0 & AF_proe == 0 & AF_pcusp == 0], unique_repIDs)
    out_d[["roemUniq"]] <- CountRepIDs(sv_summary[AF_proe > 0 & AF_pdrum == 0 & AF_pcusp == 0], unique_repIDs)
    out_d[["cuspUniq"]] <- CountRepIDs(sv_summary[AF_pcusp > 0 & AF_proe == 0 & AF_pdrum == 0], unique_repIDs)
    out_d <- as.data.table(out_d)
    out_d[, repID := unique_repIDs]

    out_fn <- str_c(out_dir, "te_sv_counts.1sp.tsv")
    cat(str_interp("Writing 1sp (i.e., unique) counts to ${out_fn}\n"))
} else if(aggr_method == "twosp" || aggr_method == "2sp"){
    # 2 species counts
    out_d <- list()
    out_d[["roem_drum"]] <- CountRepIDs(sv_summary[AF_proe > 0 & AF_pdrum > 0 & AF_pcusp == 0], unique_repIDs)
    out_d[["roem_cusp"]] <- CountRepIDs(sv_summary[AF_proe > 0 & AF_pcusp > 0 & AF_pdrum == 0], unique_repIDs)
    out_d[["drum_cusp"]] <- CountRepIDs(sv_summary[AF_pdrum > 0 & AF_pcusp > 0 & AF_proe == 0], unique_repIDs)
    out_d <- as.data.table(out_d)
    out_d[, repID := unique_repIDs]

    out_fn <- str_c(out_dir, "te_sv_counts.2sp.tsv")
    cat(str_interp("Writing 2sp counts to ${out_fn}\n"))
    
} else if(aggr_method == "threesp" || aggr_method == "3sp"){
    # 3 species counts
    out_d <- list()
    out_d[["roem_drum_cusp"]] <- CountRepIDs(sv_summary[AF_pdrum > 0 & AF_proe > 0 & AF_pcusp > 0], unique_repIDs)
    out_d <- as.data.table(out_d)
    out_d[, repID := unique_repIDs]

    out_fn <- str_c(out_dir, "te_sv_counts.3sp.tsv")
    cat(str_interp("Writing 3sp counts to ${out_fn}\n"))
} else if(aggr_method == "present" || aggr_method == "pres"){
    # "Present" counts (i.e., repeat present in an SV of a given species which may or may not be present in the other two species)
    out_d <- list()
    out_d[["drumPres"]] <- CountRepIDs(sv_summary[AF_pdrum > 0], unique_repIDs)
    out_d[["roemPres"]] <- CountRepIDs(sv_summary[AF_proe > 0] , unique_repIDs)
    out_d[["cuspPres"]] <- CountRepIDs(sv_summary[AF_pcusp > 0], unique_repIDs)
    out_d <- as.data.table(out_d)
    out_d[, repID := unique_repIDs]

    out_fn <- str_c(out_dir, "te_sv_counts.pres.tsv")
    cat(str_interp("Writing 'present' counts to ${out_fn}\n"))
} else if(aggr_method == "pressubsamp" || aggr_method == "pressub" || aggr_method == "pressamp"){
    # Present subsamp counts (i.e., subsampling "Present" SVs to match "Unique" SV numbers)
    out_d <- list()
    n <- sv_summary[AF_pdrum > 0 & AF_proe == 0 & AF_pcusp == 0, .N]
    out_d[["drumPresSubsamp"]] <- CountRepIDs(sv_summary[AF_pdrum > 0][sample.int(.N, size=n, replace=FALSE)], unique_repIDs)

    n <- sv_summary[AF_proe > 0 & AF_pdrum == 0 & AF_pcusp == 0, .N]
    out_d[["roemPresSubsamp"]] <- CountRepIDs(sv_summary[AF_proe > 0][sample.int(.N, size=n, replace=FALSE)], unique_repIDs)

    n <- sv_summary[AF_pcusp > 0 & AF_proe == 0 & AF_pdrum == 0, .N]
    out_d[["cuspPresSubsamp"]] <- CountRepIDs(sv_summary[AF_pcusp > 0][sample.int(.N, size=n, replace=FALSE)], unique_repIDs)

    out_d <- as.data.table(out_d)
    out_d[, repID := unique_repIDs]

    out_fn <- str_interp("${out_dir}/te_sv_counts.presSubsamp.${rep_i}.tsv")
    cat(str_interp("Writing subsampled 'present' counts to ${out_fn}\n"))
} else {
    stop(str_interp("Aggregation method ${aggr_method} not recognized"))
}

fwrite(out_d, file=out_fn, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
