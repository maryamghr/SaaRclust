log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source("utils/cluster_SS_reads.R")

map <- fread(snakemake@input[["aln_lib"]])
clust.to.chrom <- fread(snakemake@input[["clust_to_chrom_mapping"]])
clust.pairs <- read.table(snakemake@input[["cluster_pairs"]], header = T, stringsAsFactors = F)

fwrite(clusterSSreads(map, clust.to.chrom, clust.pairs), file = snakemake@output[[1]])


