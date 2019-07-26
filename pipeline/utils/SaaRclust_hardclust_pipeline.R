#!/usr/bin/Rscript

#This Rscript runs kmeans based hard clustering implenented in package SaaRclust.
#author: David Porubsky

args=commandArgs(TRUE)


#args <- c("aligns_k15_w-default_f0.0002_z500", "aligns_k15_w-default_f0.0002_z500/SaaRclust_results_HG00514", "100", "0.01", "30000", "TRUE", "/local/data/maryam/haploclust-css-HG00514/SaaRclust/pipeline/utils/R-packages/")
print(args)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[7]) )

suppressPackageStartupMessages(library(SaaRclust))
suppressPackageStartupMessages(library(dplyr))

output.filename <- paste0("hardClusteringResults_", args[3], "clusters.RData")

output <- runSaaRclust(inputfolder = args[1], outputfolder = args[2], num.clusters=as.numeric(args[3]), alpha = as.numeric(args[4]),  HC.only = TRUE, store.counts = FALSE, store.bestAlign = TRUE, numAlignments = as.numeric(args[5]), log.scale=args[6], outputfilename=output.filename)


