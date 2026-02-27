#!/usr/bin/env Rscript
library(data.table)

argv <- commandArgs(trailingOnly = TRUE)

bed_fn <- argv[1]
fai_fn <- argv[2]
out_fn <- argv[3]

bed <- fread(bed_fn, header = FALSE)
names(bed)[1:3] <- c("CHROM", "START", "END")
bed[, INDEX := 1:.N]
out_order <- names(bed)[]

bed[, SIZE := END - START]
min_int_size <- bed[, min(SIZE)]

fai <- fread(fai_fn,
             header = FALSE,
             select = 1:2,
             col.names = c("CHROM", "SIZE"))
fai <- fai[SIZE >= min_int_size]
chr_sizes <- fai$SIZE
names(chr_sizes) <- fai$CHROM


# Sample to get new chromosomes;
# redo those intervals whose new chromosomes
# are smaller than the interval size
m <- rep(TRUE, nrow(bed))
while (any(m)) {
    bed[m, NEW_CHROM := sample(names(chr_sizes),
                               size = .N,
                               replace = TRUE,
                               prob = unname(chr_sizes))]
    m <- bed[, chr_sizes[NEW_CHROM] < SIZE]
}

# Sample start positions
bed[, START := floor(runif(.N) * (chr_sizes[NEW_CHROM] - SIZE + 1))]
bed[, END := START + SIZE]
bed[, CHROM := NULL]
bed[, SIZE := NULL]
setnames(bed, "NEW_CHROM", "CHROM")

# Reorder columns and sort
setcolorder(bed, out_order)
setorderv(bed, c("CHROM", "START"))

fwrite(bed,
       file = out_fn,
       sep = "\t",
       col.names = FALSE,
       row.names = FALSE,
       quote = FALSE)
