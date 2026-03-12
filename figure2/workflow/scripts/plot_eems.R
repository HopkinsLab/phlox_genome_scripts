#!/usr/bin/env Rscript

library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(rworldxtra)

argv <- commandArgs(trailingOnly=TRUE)

mcmcpath = argv[1]
plotpath = argv[2]

projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"

eems.plots(mcmcpath, plotpath, 
           longlat = FALSE,
           add.grid = TRUE,
           add.demes = TRUE,
           projection.in = projection_none,
           projection.out = projection_mercator,
           add.map = TRUE,
           col.map = "black")

