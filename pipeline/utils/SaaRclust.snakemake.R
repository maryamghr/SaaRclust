log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

.libPaths(c(.libPaths(), 'utils/R-packages/'))
library(dplyr)
library(SaaRclust)
library(data.table)
library(Rsamtools)
library(doParallel)
library(matrixStats)
library(assertthat)
library(biovizBase)
source('../R/timedMessage.R')
source('../R/import.R')
source('../R/hardClust.R')
source('../R/SaaRclust.R')
source('../R/wrapper.R')
source('../R/EMclust.R')
source('../R/calcProbs.R')
source('../R/helperFuctions.R')
source('../R/dumped_functions.R')

inputfolder <- dirname(snakemake@input[["bam"]][1])
outputfolder <- dirname(dirname(snakemake@output[["hard_clust"]]))
input_type <- snakemake@params[["input_type"]]
input.alignment.files <- snakemake@input[["bam"]]
ref.aln.bam <- snakemake@input[["unitigs_bam"]]
num.clusters <- as.numeric(snakemake@params[["num_clusters"]])
EM.iter <- as.numeric(snakemake@params[["EMiter"]])
numAlignments <- as.numeric(snakemake@params[["num_alignments"]])
hardclust.file <- snakemake@output[["hard_clust"]]
softclust.file <- snakemake@output[["soft_clust"]]
MLclust.file <- snakemake@output[["ML_clust"]]
ss.clust.file <- snakemake@output[["ss_clust"]]
clust.pairs.file <- snakemake@output[["clust_pairs"]]
wc.cells.file <- snakemake@output[["wc_cells_clusters"]]
numCPU <- snakemake@threads[[1]]

clust <- runSaaRclust(inputfolder=inputfolder, outputfolder=outputfolder, 
                      num.clusters=num.clusters, EM.iter=EM.iter, 
                      numAlignments=numAlignments, hardclust.file=hardclust.file, 
                      softclust.file=softclust.file, MLclust.file=MLclust.file, 
                      ss.clust.file=ss.clust.file, clust.pairs.file=clust.pairs.file, 
                      wc.cells.file=wc.cells.file, ref.aln.bam=ref.aln.bam, numCPU=numCPU)
