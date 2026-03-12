#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)


argv <- commandArgs(trailingOnly = TRUE)
in_prefix <- argv[1]
k_vals <- argv[2]
met_fn <- argv[3]
out_prefix <- argv[4]

# in_prefix <- "merged.roem.chrOnly_1.minGQ_30.admix"
# k_vals <- "1,2,3,4,5,6,7,8"
# met_fn <- "metadata.tsv"
# out_prefix <- "foo"

k_vals <- as.integer(str_split(k_vals, ",")[[1]])
k_vals <- k_vals[!is.na(k_vals)]

# Read in metadata
met_d <- fread(met_fn, header=TRUE, sep="\t", na.strings=c("NA", "#N/A"))

# Read in CV report
cv_d <- fread(str_interp("${in_prefix}.CV_report.txt"), sep=":", header=FALSE,
              col.names=c("K", "CV_error"))
cv_d[, K := as.integer(str_replace(str_replace(K, "CV error \\(K=", ""), "\\)", ""))]

# Get minimum k value
min_k <- cv_d[which.min(CV_error), K]

# Read in admixture data
samps <- fread(str_interp("${in_prefix}.fam"), header=FALSE)$V2
admix_d <- list()
for(k in k_vals){
  fn <- str_interp("${in_prefix}.${k}.Q")
  admix_d[[k]] <- fread(fn, header=FALSE)
  admix_d[[k]][, IID := samps]
  admix_d[[k]] <- melt(admix_d[[k]], id.vars = "IID", variable.name="pop_infer")
  admix_d[[k]][, K := k]
  if(k == min_k){
    admix_d[[k]][, K_label := str_interp("* K = ${k} *")] # Star if minimum
  } else {
    admix_d[[k]][, K_label := str_interp("K = ${k}")]
  }
}
admix_d <- rbindlist(admix_d)
admix_d[, pop_infer := str_replace(pop_infer, "V", "p")]

admix_long <- merge(admix_d, met_d, by.x="IID", by.y="Sample_ID")
admix_long <- admix_long[!is.na(species) & !is.na(region)] # Remove unknown species
admix_long[, species_region := str_c(species, region, sep="_")]

# Reorder species based on species + region
samp_ord <- met_d[Sample_ID %in% unique(admix_long$IID), 
               .(Sample_ID, species_region = str_c(species, region, sep="_"))]
setorder(samp_ord, species_region)
admix_long[, IID := factor(IID, levels=samp_ord$Sample_ID)]


# Get vertical line positions
x_annot <- which(!duplicated(samp_ord$species_region))[-1] - 0.5
# Get text annotations

# Order K labels
k_ord <- admix_long[, .(K_label = K_label[1]), by=K]
setorder(k_ord, K)
k_levels <- k_ord$K_label
names(k_levels) <- as.character(k_ord$K)
admix_long[, K_label := factor(K_label, levels=unname(k_levels[as.character(k_vals)]))]

#######################
# Admixture plots
plt_d <- admix_long[K != 1,]

x_breaks <- levels(plt_d$IID)[c(1, x_annot+0.5)]

x_ticklabs <- samp_ord[!duplicated(species_region), species_region]
x_ticklabs <- str_replace(x_ticklabs, "_", "\n")


p <- ggplot(plt_d, aes(x=IID, y=value, fill=pop_infer)) + 
  geom_bar(stat="identity", width=1) +
  facet_wrap(vars(K_label), nrow=length(k_vals)) + theme_bw() +
  ylab("") +
  scale_x_discrete(name="", breaks=x_breaks, labels=x_ticklabs) +
  theme(axis.text.x = element_text(size=7, hjust=0),
        axis.ticks.x = element_blank())
for(x in x_annot){
  p <- p + geom_vline(xintercept = x, color="black")
}
p <- set_palette(p, "d3")
# print(p)

########################################################
# Save admixture plot
GetPlotDimensions <- function(plt_d, px_per_k=118.5, px_per_samp=2.5, w.pad=160, h.pad=50){
  h <- length(unique(plt_d$K))   * px_per_k    + h.pad
  w <- length(unique(plt_d$IID)) * px_per_samp + w.pad
  return(c("height"=round(h), "width"=round(w)))
}

p_dim <- GetPlotDimensions(plt_d)
out_fn <- str_interp("${out_prefix}.oBy_species_region.png")
ggsave(p, file=out_fn, width=p_dim["width"], height=p_dim["height"], units="px", dpi=96)


########################################################
# CV error plot
cv_d[, IS_MIN := FALSE]
cv_d[which.min(CV_error), IS_MIN := TRUE]
p <- ggplot(cv_d, aes(x=K, y=CV_error)) + 
  geom_line(color="#444444") + geom_point(aes(color=IS_MIN)) + 
  theme_pubr() + ylab("CV error") + 
  scale_color_manual(values=c("TRUE"="coral2", "FALSE"="#444444")) +
  guides(color="none")
# print(p)

out_fn <- str_interp("${out_prefix}.CV_error.png")
ggsave(p, file=out_fn, width=5, height=4, units="in")

