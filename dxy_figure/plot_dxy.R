#!/usr/bin/env Rscript
library(data.table)
library(stringr)
library(ggpubr)
library(patchwork)

source("../utils/genome_plots.R")

dat <- list()
dat[["WGS"]] <- fread("merged.SUMMED.mpileup.DP_10.all_dxy.txt", header = TRUE)
dat[["RRS"]] <- fread("multi.cusp_drum_roem_pilo.RRS_dxy.txt", header = TRUE)

#####

plt_d <- dat$RRS[pop2 == "pilo", .(DXY = sum(count_diffs) / sum(count_comparisons)), by = .(pop1, pop2, MIN_SAMP, MIND)]
print(plt_d[MIN_SAMP == 0.7][MIND == 0.25])

plt_d[, GROUP := str_c(pop1, pop2, sep = "-")]
plt_d[, MIN_SAMP := as.numeric(MIN_SAMP)]
plt_d[, FACET := str_c("Max per-sample missingness: ", MIND), by = seq(.N)]
use_colors <- Load4SpPal()

setnames(plt_d, "pop1", "Species")
p1 <- ggplot(plt_d, aes(x = MIN_SAMP, y = DXY, color = Species, group = Species)) +
    geom_line() +
    geom_point() +
    xlab("Minimum proportion of samples\nto process locus") +
    ylab("dxy vs. P. pilosa") +
    # geom_point(shape = 21) +
    # scale_fill_manual(values = use_colors) +
    scale_color_manual(values = use_colors) +
    facet_wrap(vars(FACET), ncol = 1) +
    ggtitle("A") +
    theme_pubr() +
    theme(legend.position = "none")
print(p1)

#####
plt_d <- dat$WGS[pop2 == "pilo", .(Y = median(avg_dxy)), by = .(Species = pop1)]

use_colors <- Load4SpPal()

p2 <- ggbarplot(plt_d, x = "Species", y = "Y", fill = "Species") +
    scale_fill_manual(values = use_colors) +
    xlab("Species") +
    ylab("dxy vs. P. pilosa") +
    ggtitle("B") +
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
