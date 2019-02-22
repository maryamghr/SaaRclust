log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source("utils/mappingClusterNames.R")

soft.clust.probs <- snakemake@input[[1]] #"soft_probs/NA12878_WashU_PBreads_soft_probs.data.gz"
cluster.mapping <- snakemake@output[[1]] #"soft_probs/cluster_mapping.data"

clusterMapping(soft.clust.probs, cluster.mapping)
