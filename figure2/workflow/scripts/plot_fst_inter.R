#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)
library(ggpubr)

# Arguments
argv <- commandArgs(trailingOnly = TRUE)
# argv <- c("test.list", "1000000", "results/tmp")

fn_list <- fread(argv[1], header=FALSE)$V1
wsize <- as.numeric(argv[2])
out_dir <- argv[3]
s <- "fst"
yvar <- "avg_wc_fst"

# Functions
ParseFilenameForFilters <- function(fn, filter_prefixes=c("min_samples_per_pop~", "mind_")){
    split_path <- str_split(fn, "/")[[1]]
    out_strs <- c()
    for(fp in filter_prefixes){
        new_str <- split_path[str_starts(split_path, fp)]
        # new_str <- str_c(str_remove(fp, ".$"), str_remove(tmp, fp), sep=": ")
        out_strs <- c(out_strs, new_str)
    }
    return(out_strs)
}

AbbrevBaseLen <- function(x){
    if(x < 1){
        stop("base length 'x' must be positive integer")
    }
    pwr <- floor(log10(x))

    if(pwr < 3){
        out_str <- str_interp("$[d]{x} bp")
    } else if(pwr < 6){
        out_str <- str_interp("$[0.1f]{x / (10^3)} kb")
    } else if(pwr < 9){
        out_str <- str_interp("$[0.1f]{x / (10^6)} Mb")
    } else if(pwr < 12){
        out_str <- str_interp("$[0.1f]{x / (10^9)} Gb")
    } else {
        stop("base length 'x' is too large to abbreviate")
    }
    return(out_str)
}

########################################
# Read in data
dat <- list()
for(fn in fn_list){ # Iterate over files in filelist
    dat[[fn]] <- fread(fn, header=TRUE)[!is.na(get(yvar))] # Remove NAs
    filt_strs <- ParseFilenameForFilters(fn)
    dat[[fn]][, FILENAME := fn]
    dat[[fn]][, FILTER1 := filt_strs[1]]
    dat[[fn]][, FILTER2 := filt_strs[2]]
}
dat <- rbindlist(dat)
dat[, POS := window_pos_2 / 1e6 ] # Add


########################################
# Plot along genome
cvar <- "Comparison"
use_palette <- "lancet"
use_ylab <- "Average Fst"

dat[, Comparison := str_c(pop1, pop2, sep=" x ")]
xbreaks <- seq(0, dat[, max(POS)], by=50)
for(fn in fn_list){
    filt_strs <- ParseFilenameForFilters(fn)
    prefix <- str_c(out_dir, "/", str_remove(basename(fn), "\\.txt$"), 
                    ".", filt_strs[1], ".", filt_strs[2])
        
    p <- ggline(dat[FILENAME == fn], x="POS", y=yvar, plot_type="l", numeric.x.axis=TRUE, 
                size=0.2, color=cvar, palette=use_palette) + 
        facet_wrap(vars(chromosome), nrow=7) +
        scale_x_continuous(name="Coordinate (Mb)", breaks=xbreaks) +
        ylab(use_ylab)
    out_fn <- str_interp("${prefix}.byChr.pdf")
    ggsave(p, filename=out_fn, width=8, height=12)
}

########################################
# Plot genome-wide pi
xvar <- "Comparison"
gnm_d <- dat[, .(Y = weighted.mean(avg_wc_fst, no_snps)), by=c(xvar, "FILTER1", "FILTER2")]

nr <- length(unique(gnm_d[,FILTER1]))
nc <- length(unique(gnm_d[,FILTER2]))
# Genome values (barplot)
p <- ggbarplot(gnm_d, x=xvar, y="Y", fill=xvar, palette=use_palette) +
     ylab(str_interp("Genome-wide ${s}")) +
     facet_grid(FILTER1 ~ FILTER2)
out_fn <- str_interp("${out_dir}/${s}.genome_wide.bar.pdf")
ggsave(p, filename=out_fn, width=3*nc, height=3*nr)

# Quantile values (boxplot)
p <- ggboxplot(dat, x=xvar, y=yvar, fill=xvar, palette=use_palette, alpha=0.6) +
     ylab(str_interp("${s} per ${AbbrevBaseLen(wsize)}")) +
     facet_grid(FILTER1 ~ FILTER2)
out_fn <- str_interp("${out_dir}/${s}.genome_wide.box.pdf")
ggsave(p, filename=out_fn, width=3*nc, height=3*nr)
