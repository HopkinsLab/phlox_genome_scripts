#!/usr/bin/env Rscript
library(rmarkdown)

argv <- commandArgs(trailingOnly=TRUE)

render(argv[1])
