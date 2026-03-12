#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)
library(ggpubr)
library(geosphere)

# Arguments
argv <- commandArgs(trailingOnly = TRUE)
# argv <- c("test.list", "1000000", "metadata.tsv", "results/tmp")

fn_list <- fread(argv[1], header=FALSE)$V1
wsize <- as.numeric(argv[2])
met_fn <- argv[3]
out_dir <- argv[4]
s <- "fst"

################ 
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

ParseFilenameForSpecies <- function(fn){
    split_path <- str_split(fn, "/")[[1]]
    for(sp in c("cusp", "drum", "roem")){
        if(sp %in% split_path){
            return(sp)
        }
    }
    stop(str_interp("${fn} does not contain a recognizeable species name."))
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

met_d <- fread(met_fn, header=TRUE, sep="\t")

dat <- list()
for(fn in fn_list){ # Iterate over files in filelist
    dat[[fn]] <- fread(fn, header=TRUE)[!is.na(avg_wc_fst), .(FST = weighted.mean(avg_wc_fst, w=no_snps)), by=.(pop1, pop2)] # Remove NAs
    filt_strs <- ParseFilenameForFilters(fn)
    sp <- ParseFilenameForSpecies(fn)
    dat[[fn]][, MD5 := digest::digest(filt_strs, algo="md5")]
    dat[[fn]][, SPECIES := sp]
    dat[[fn]][, FILTER1 := filt_strs[1]]
    dat[[fn]][, FILTER2 := filt_strs[2]]
}
dat <- rbindlist(dat)

########################################
# Compute km distance between populations

compute_haversine <- function(p1, p2, d){
    coord1 <- d[(p1), c(long[1], lat[1])]
    coord2 <- d[(p2), c(long[1], lat[1])]
    return(distHaversine(coord1, coord2))
}

super_pops <- dat[, unique(c(pop1, pop2))]

pop_d <-  list()
pop_d[[1]] <- met_d[, .(lat=lat[1], long=long[1]), by=Population]
setkey(pop_d[[1]], Population)
for(x in setdiff(super_pops, pop_d[[1]]$Population)){
    y <- str_split(x, "_")[[1]][1]
    pop_d[[x]] <- pop_d[[1]][(y),][1]
    pop_d[[x]][, Population := x]
}
pop_d <- rbindlist(pop_d)
setkey(pop_d, Population)

dist_d <- dat[, .(INDEX=1), by=.(pop1,pop2)]
dist_d[, INDEX := .I]
dist_d[, DIST := compute_haversine(pop1, pop2, pop_d) / 1e3, by=INDEX]
dat <- merge(dat, dist_d[, .(pop1, pop2, DIST)], by=c("pop1", "pop2"))

########################################
# Linear model

OLS <- function(d){
    m <- lm(FST ~ DIST, d)
    pval <- coef(summary(m))[2,4]
    b0 <- coef(m)[1]
    b1 <- coef(m)[2]

    pred_d <- d[, .(DIST = seq(min(DIST), max(DIST), length.out=100))]
    tmp <- predict(m, newdata=pred_d, interval="confidence")
    pred_d[, FST := tmp[, "fit"]]
    pred_d[, LWR := tmp[, "lwr"]]
    pred_d[, UPR := tmp[, "upr"]]
    return(list(pred=pred_d, pval=pval, b0=b0, b1=b1))
}

line_d <- dat[, OLS(.SD)[["pred"]], by=.(SPECIES, FILTER1, FILTER2)]
annot_d <- dat[, OLS(.SD)[c("pval", "b0", "b1")], by=.(SPECIES, FILTER1, FILTER2)]
annot_d[, TEXT := sprintf("\n p = %0.1e\n b0 = %0.1fe-3\n b1 = %0.3fe-3", pval, b0*1e3, b1*1e3)]

########################################
# Plot FST vs DIST
pal <- get_palette("d3", 3)
names(pal) <- dat[, sort(unique(SPECIES))]

nr <- length(unique(dat[,FILTER1]))
nc <- length(unique(dat[,FILTER2]))

for(species in names(pal)){
    p <- ggscatter(dat[SPECIES == species], x="DIST", y="FST", color=pal[species], alpha=0.3) +
        geom_ribbon(aes(x=DIST, ymin=LWR, ymax=UPR), data=line_d[SPECIES == species], fill="gray", color=NA, alpha=0.5) +
        geom_line(aes(x=DIST, y=FST), data=line_d[SPECIES == species], color=pal[species]) +
        geom_text(aes(x=-Inf, y=Inf, hjust=0, vjust=1, label=TEXT), data=annot_d[SPECIES == species]) +
        xlab("Distance (km)") +
        ylab(str_interp("Average ${s}")) +
        facet_grid(FILTER1 ~ FILTER2) + 
        theme_bw()
    out_fn <- str_interp("${out_dir}/${s}.${species}.scatter.vs_dist.pdf")
    # print(out_fn)
    ggsave(p, filename=out_fn, width=3*nc, height=3*nr)
}
