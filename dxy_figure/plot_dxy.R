#!/usr/bin/env Rscript
library(data.table)
library(stringr)
library(ggpubr)
library(patchwork)

source("../utils/genome_plots.R")

#### Load data
dat <- list()
dat[["RRS"]] <- list()
for (minsamp in c("0.1", "0.5", "0.7")) {
    for (mind in c("0.25", "0.5", "0.75", "0.9")) {
        x <- str_c(minsamp, "::", mind)
        fn <- str_interp("RRS_data/multi.cusp_drum_roem_pilo.minsamp_${minsamp}.mind_${mind}.RRS_dxy.txt")
        dat[["RRS"]][[x]] <- fread(fn, header = TRUE)
        dat[["RRS"]][[x]][, MIN_SAMP := as.numeric(minsamp)]
        dat[["RRS"]][[x]][, MIND     := as.numeric(mind)]
    }
}
dat[["RRS"]] <- rbindlist(dat[["RRS"]])

dat[["WGS"]] <- list()
for (i in seq_len(5)) {
    fn <- str_interp("WGS_data/merged.CONCAT.mpileup.rep${i}.DP_10.all_dxy.txt")
    dat[["WGS"]][[i]] <- fread(fn, header = TRUE)
    dat[["WGS"]][[i]][, REP := i]
}
dat[["WGS"]] <- rbindlist(dat[["WGS"]])

##### Filter for just chromosomes
dat$RRS <- dat$RRS[no_sites > 0][str_starts(chromosome, "chr")]
dat$WGS <- dat$WGS[no_sites > 0]

##### Summarize data
sum_dat <- list()
sum_dat[["WGS"]] <- dat$WGS[
    ,
    .(no_sites          = sum(as.numeric(no_sites)),
      count_diffs       = sum(as.numeric(count_diffs)),
      count_comparisons = sum(as.numeric(count_comparisons)),
      count_missing     = sum(as.numeric(count_missing))),
    by = .(pop1, pop2, REP)
]
sum_dat$WGS[, avg_dxy := count_diffs / count_comparisons]
setcolorder(sum_dat$WGS, neworder = intersect(names(dat$WGS), names(sum_dat$WGS)))

sum_dat[["RRS"]] <- dat$RRS[
    ,
    .(no_sites          = sum(as.numeric(no_sites)),
      count_diffs       = sum(as.numeric(count_diffs)),
      count_comparisons = sum(as.numeric(count_comparisons)),
      count_missing     = sum(as.numeric(count_missing))),
    by = .(pop1, pop2, MIN_SAMP, MIND)
]
sum_dat$RRS[, avg_dxy := count_diffs / count_comparisons]
setcolorder(sum_dat$RRS, neworder = intersect(names(dat$RRS), names(sum_dat$RRS)))

##### Resample data (for CIs)
BootstrapWindows <- function(numer, denom, nboot = 1000) {
    dxy <- colSums(replicate(sample(numer, replace = TRUE), n = nboot)) /
        colSums(replicate(sample(denom, replace = TRUE), n = nboot))

    res <- as.list(quantile(dxy, probs = c(0.025, 0.5, 0.975)))
    names(res) <- c("lwr", "med", "upr")
    res
}

ci_dat <- list()
set.seed(4558)

ci_dat[["WGS"]] <- dat$WGS[
    ,
    BootstrapWindows(count_diffs, count_comparisons),
    by = .(pop1, pop2, REP)
]

ci_dat[["RRS"]] <- dat$RRS[
    ,
    BootstrapWindows(count_diffs, count_comparisons),
    by = .(pop1, pop2, MIN_SAMP, MIND)
]

##### Plot

plt_d <- merge(
    sum_dat$RRS[pop2 == "pilo"],
    ci_dat$RRS[pop2 == "pilo"],
    by = c("pop1", "pop2", "MIN_SAMP", "MIND")
)

print(plt_d[MIN_SAMP == 0.7][MIND == 0.25])

setnames(plt_d, "pop1", "Species")
plt_d[, Species := factor(Species, levels = c("roem", "drum", "cusp"))]
plt_d[, OFFSET := (2 - as.integer(Species)) * 0.15]
plt_d[, MIN_SAMP := factor(as.character(MIN_SAMP))]
plt_d[, FACET := str_c("Max per-sample missingness: ", MIND), by = seq(.N)]
use_colors <- Load4SpPal()

p1 <- ggplot(
    plt_d,
    aes(x = as.integer(MIN_SAMP) + OFFSET,
        y = avg_dxy,
        ymin = lwr,
        ymax = upr,
        color = Species,
        group = Species)) +
    geom_point() +
    geom_linerange() +
    scale_x_continuous(
        name = "Minimum proportion of samples\nto process locus",
        breaks = seq_len(3),
        labels = levels(plt_d$MIN_SAMP)
    ) +
    ylab("dxy vs. P. pilosa") +
    scale_color_manual(values = use_colors) +
    facet_wrap(vars(FACET), ncol = 1) +
    ggtitle("A") +
    theme_pubr() +
    theme(legend.position = "none")
print(p1)

#####
plt_d <- merge(
    sum_dat$WGS[pop2 == "pilo", .(Y = median(avg_dxy)), by = .(Species = pop1)],
    ci_dat$WGS[pop2 == "pilo", .(LWR = median(lwr), UPR = median(upr)), by = .(Species = pop1)],
    by = "Species"
)

use_colors <- Load4SpPal()

p2 <- ggplot(plt_d, aes(x = Species, y = Y, ymin = LWR, ymax = UPR, color = Species)) +
    geom_point() +
    geom_linerange() +
    scale_color_manual(values = use_colors) +
    xlab("Species") +
    ylab("dxy vs. P. pilosa") +
    ggtitle("B") +
    theme_pubr() +
    theme(legend.position = "right")
print(p2)

print(p1 + p2)

fn <- "Supplement.dxy.pdf"
ggsave(
    p1 + p2,
    file = fn,
    width = 9,
    height = 6,
    device = cairo_pdf,
    family = "Arial"
)
