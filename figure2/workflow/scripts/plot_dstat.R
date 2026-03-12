#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)
library(ggpubr)


argv <- commandArgs(trailingOnly=TRUE)
# argv <- c(  "/n/holyscratch01/hopkins_lab/Users/felixw/demog_annuals/ddRAD_proc/align_pilosa/results/multi/cusp_drum_roem_pilo/res_stacks/min_samples_per_pop~0.1/min_samples_overall~0.1/mind_0.25/dstat/populations.snps.lFilt.iFilt.patt_counts.tsv",
#             "metadata.tsv",
#             "results/tmp/dstat_test")


fn <- argv[1]
out_prefix <- argv[2]
# met_fn <- argv[2]


###### Load data

dat <- fread(fn, header=TRUE, sep="\t")
sp_names <- str_subset(names(dat), "0", negate=TRUE)
dat[, D := (`0110` - `1010`) / (`0110` + `1010`)]

# met <- fread(met_fn, header=TRUE, sep="\t")


###### Histogram
p <- gghistogram(dat, x="D", bins=nclass.FD(dat$D), fill="dodgerblue")

out_fn <- str_interp("${out_prefix}.hist.pdf")
cat(out_fn, "\n", sep="")
ggsave(p, filename=out_fn, width=4, height=3)

###### Correlation between "(BBAA - max(ABBA, BABA)) / BBAA" and "D"
dat[, BBAA_diff := (`1100` - max(`0110`, `1010`)) / `1100`, by=1:nrow(dat)]
dat[, pos_D := D > 0]
p <- ggscatter( dat, x="BBAA_diff", y="D", color="pos_D", palette="Dark2",
                xlab="(BBAA - max(ABBA, BABA)) / BBAA", ylab="D") +
        geom_vline(xintercept=0) +
        geom_hline(yintercept=0)

out_fn <- str_interp("${out_prefix}.BBAA_diff.scatter.pdf")
cat(out_fn, "\n", sep="")
ggsave(p, filename=out_fn, width=4, height=4)


# ###### Where are negative and positive drum values in Texas?
# tmp <- list()
# texas_d <- list()
# for(sp in sp_names){
#     tmp[[sp]] <- dat[, .(Sample_ID=unique(get(sp)), Species=sp), by=.(pos_D = D > 0)]
#     texas_d[[sp]] <- as.data.table(map_data("state", "texas"))
#     texas_d[[sp]][, Species := sp]
# }
# tmp <- rbindlist(tmp)
# texas_d <- rbindlist(texas_d)

# plt_d <- merge(tmp, met, by="Sample_ID")[, .(Sample_ID, lat, long, pos_D, Species)]

# use_col <- scales::viridis_pal()(2)
# names(use_col) <- c("TRUE", "FALSE")

# p <- ggplot() +
#         geom_polygon(aes(long, lat, group = group), data=texas_d, fill = "grey90", colour = "grey30") +
#         geom_point(aes(long, lat, fill=pos_D), data=plt_d, color="grey60", shape=21, alpha=0.5, position=position_jitter(width=0.1)) +
#         scale_fill_manual(values=use_col) + 
#         facet_wrap(vars(Species), nrow=1) +
#         coord_quickmap(xlim=c(-100, -95), ylim=c(28, 33)) + 
#         theme_bw() + theme(legend.position="top")

# out_fn <- str_interp("${out_prefix}.map.pdf")
# cat(out_fn, "\n", sep="")
# ggsave(p, filename=out_fn, width=12, height=4)
