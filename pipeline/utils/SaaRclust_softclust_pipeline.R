#!/usr/bin/Rscript

#This Rscript runs EM (Expectation Maximization) based soft clustering algorithm implenented in package SaaRclust.
#author: David Porubsky

args=commandArgs(TRUE)
print(args)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[11]) )

suppressPackageStartupMessages(library(SaaRclust))

output <- SaaRclust(minimap.file=args[1], outputfolder=args[2], num.clusters=args[3], EM.iter=args[4], alpha=as.numeric(args[5]), minLib=as.numeric(args[6]), upperQ=as.numeric(args[7]), logL.th=args[8], theta.constrain=FALSE, HC.input=args[9], log.scale=args[10])

# testing
args <- c("aligns_k15_w1_f0.05_z8/NA12878_WashU_PBreads_chunk000.maf.gz", "aligns_k15_w1_f0.05_z8/SaaRclust_results_NA12878_WashU_PBreads", "47", "2", "0.01", "1", "0.95", "1", "aligns_k15_w1_f0.05_z8/SaaRclust_results_NA12878_WashU_PBreads/Clusters/hardClusteringResults.RData", "TRUE", "/MMCI/TM/scratch/maryam/clustering/SaaRclust/pipeline/utils/R-packages/")
