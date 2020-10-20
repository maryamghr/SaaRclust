#!/usr/bin/Rscript

#This Rscript runs kmeans based hard clustering implenented in package SaaRclust.
#author: David Porubsky, Maryam Ghareghani

args=commandArgs(TRUE)
print(version)

.libPaths( c( .libPaths(), args[13]) )

library(data.table)
suppressPackageStartupMessages(library(SaaRclust))
suppressPackageStartupMessages(library(dplyr))
library(Rsamtools)
library(doParallel)
library(matrixStats)
source('../R/timedMessage.R')

clust <- runSaaRclust(inputfolder = args[1], outputfolder = args[2], input_type=args[3], num.clusters=as.numeric(args[4]),
                      EM.iter=args[8], alpha=as.numeric(args[5]), minLib=as.numeric(args[9]), upperQ=as.numeric(args[10]),
                      logL.th=args[11], theta.constrain=FALSE, store.counts=FALSE, store.bestAlign=TRUE, numAlignments=as.numeric(args[6]),
                      HC.only=FALSE, log.scale=args[7], numCPU=as.numeric(args[12]))
