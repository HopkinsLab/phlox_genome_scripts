#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)
library(ggpubr)

argv <- commandArgs(trailingOnly=TRUE)
# in_fn <- "results/multi/cusp_drum_roem_pilo/res_stacks/min_samples_per_pop~0.7/min_samples_overall~0.7/mind_0.75/freq/populations.all.lFilt.iFilt.freq.tsv"
# out_prefix <- "test"
in_fn <- argv[1]
out_prefix <- argv[2]

dat <- fread(cmd=str_interp("grep -v , ${in_fn}"), sep="\t", header=TRUE)
dat[, INDEX := .I]

all_sp <- c("cusp", "drum", "roem", "pilo")

# Calculate AF
for(sp in all_sp){
    an <- str_c("AN_", sp)
    ac <- str_c("AC_", sp)
    new_col <- str_c("AF_", sp)
    dat[, eval(new_col) := get(ac) / get(an)]
}

# Calculate subsampled AF
CalcAF <- function(an, ac_hom, ac_het, gt_sampsize=NULL){
  if(is.null(gt_sampsize) || (gt_sampsize == an/2)){
    af <- (ac_hom + ac_het) / an
  } else {
    af <- sum(sample(c(rep(1, ac_het), rep(2, ac_hom/2), rep(0, an/2 - ac_het - ac_hom/2)), size=gt_sampsize, replace=FALSE)) / (2*gt_sampsize)
  }
  return(af)
}

dat[, MIN_CALLED := min(AN_cusp, AN_drum, AN_roem)/2, by=INDEX]
min_sampsize_range <- range(dat$MIN_CALLED)

for(sp in all_sp[1:3]){
    an <- str_c("AN_", sp)
    ac_hom <- str_c("AC_Hom_", sp)
    ac_het <- str_c("AC_Het_", sp)
    new_col <- str_c("AF_Subsamp_", sp)
    dat[, eval(new_col) := CalcAF(get(an), get(ac_hom), get(ac_het), gt_sampsize=MIN_CALLED), by=INDEX]
}
dat[, AF_Subsamp_pilo := AF_pilo]



# Plot
nbins <- 20
binw <- 1/nbins
x_ticks <- seq(0, 1, by=binw)
y_ticks <- 10**c(1:5)
use_colors <- get_palette("d3", 4)
names(use_colors) <- all_sp

for(pre in c("AF", "AF_Subsamp")){
    plt_d <- melt(  dat, measure.vars=str_c(pre, "_", all_sp),
                    variable.name="SPECIES", value.name="AF")
    plt_d[, SPECIES := str_remove(SPECIES, str_c(pre, "_"))]
    
    scientific_10 <- function(x) {
        parse(text=gsub("1e\\+", "10^", scales::scientific_format()(x)))
    }
    

    if(pre == "AF"){
        use_title <- str_interp("${dat[,.N]} SNPs")
    } else {
        use_title <- str_interp("${dat[,.N]} SNPs\nAnnuals' GTs subsampled to smallest\nspecies sample size [${min_sampsize_range[1]}, ${min_sampsize_range[2]}]")
    }

    p <- gghistogram(plt_d, x="AF", fill="SPECIES", alpha=0.9, binwidth=binw, color=NA, 
                    position=position_dodge2(width=1/nbins, preserve="single", padding=0), 
                    xlab="Allele Frequency",
                    title=use_title) +
        scale_fill_manual(values=use_colors) +
        geom_vline(xintercept=seq(binw, 1, by=binw) - (binw/2), color="gray", linetype=2) +
        scale_x_continuous(breaks=x_ticks, labels=sprintf("%0.2f", x_ticks)) +
        scale_y_continuous(breaks=y_ticks, labels=scientific_10, trans=scales::pseudo_log_trans(base=10)) +
        theme(axis.text.x = element_text(angle=45, hjust=1))

    out_fn <- str_interp("${out_prefix}.${pre}.hist.pdf")
    cat(out_fn, "\n", sep="")
    ggsave(p, file=out_fn, width=6, height=4)

}
