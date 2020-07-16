#!/usr/bin/Rscript

#This Rscript runs EM (Expectation Maximization) based soft clustering algorithm implenented in package SaaRclust.
#author: David Porubsky

source('../R/timedMessage.R')

args=commandArgs(TRUE)
print(args)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[13]) )

suppressPackageStartupMessages(library(SaaRclust))
library(dplyr)
library(data.table)

print(SaaRclust)

source('../R/SaaRclust.R')

print(SaaRclust)

output <- SaaRclust(minimap.file=args[1], outputfolder=args[2], num.clusters=args[3], EM.iter=args[4], alpha=as.numeric(args[5]), minLib=as.numeric(args[6]), upperQ=as.numeric(args[7]), logL.th=args[8], theta.constrain=FALSE, HC.input=args[9], log.scale=args[10], filter.soft.clust.input=args[11], filter.ss.file=args[12])

# testing
args <- c("aligns_k15_w1_f0.1_z500/HG00733_chunk028.maf.gz", "aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733", 47, 2, 0.01, 1, 1, 1,           "aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/hardClusteringResults_80clusters.RData", "TRUE", "FALSE", "$(pwd)/utils/R-packages/")
minimap.file=args[1]
outputfolder=args[2]
num.clusters=args[3]
EM.iter=as.numeric(args[4])
alpha=as.numeric(args[5])
minLib=as.numeric(args[6])
upperQ=as.numeric(args[7])
logL.th=args[8]
theta.constrain=FALSE
HC.input=args[9]
log.scale=args[10]
filter.soft.clust.input=args[11]
ss.names.filter=args[12]
theta.param=NULL
pi.param=NULL
cellNum=NULL
store.counts=TRUE

