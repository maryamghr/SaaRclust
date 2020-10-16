#!/usr/bin/Rscript

#This Rscript runs kmeans based hard clustering implenented in package SaaRclust.
#author: David Porubsky

args=commandArgs(TRUE)
print(version)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")

#args <- c("aligns_k15_w-default_f0.0002_z500", "aligns_k15_w-default_f0.0002_z500/SaaRclust_results_HG00514", "100", "0.01", "30000", "TRUE", "/local/data/maryam/haploclust-css-HG00514/SaaRclust/pipeline/utils/R-packages/")
print(args)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[7]) )

print('library(data.table)')
library(data.table)
print('library(SaaRclust)')
suppressPackageStartupMessages(library(SaaRclust))
print(filterInput)
print('library(dplyr)')
suppressPackageStartupMessages(library(dplyr))
print('running hard clustering')
output.filename <- paste0("hardClusteringResults_", args[3], "clusters.RData")

#hard.clust <- runSaaRclust(inputfolder = args[1], outputfolder = args[2], input_type=args[3], num.clusters=as.numeric(args[4]), 
#                           alpha = as.numeric(args[5]),  HC.only = TRUE, store.counts = FALSE, store.bestAlign = TRUE, numAlignments = as.numeric(args[6]), 
#                           log.scale=args[7], outputfilename=output.filename)

hard.clust <- runSaaRclust(inputfolder = args[1], outputfolder = args[2], input_type=args[3], num.clusters=as.numeric(args[4]), 
                         EM.iter=args[8], alpha=as.numeric(args[5]), minLib=as.numeric(args[9]), upperQ=as.numeric(args[10]), 
                         logL.th=args[11], theta.constrain=FALSE, store.counts=FALSE, store.bestAlign=TRUE, numAlignments=as.numeric(args[6]), 
                         HC.only=FALSE, log.scale=args[7], outputfilename=output.filename)
  