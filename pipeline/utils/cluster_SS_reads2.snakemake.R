log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source("utils/clustering_functions.R")

cov <- lapply(snakemake@input[["ss_clust_count"]], fread)
clust.partners <- fread(snakemake@input[["clust_partners"]])

d <- mergeAllCoverages(cov, clust.partners)
d <- callMajorityVoteCluster(d)
evaluateSSclustering(d)

fwrite(d[, .(clust.forward, name)], file=snakemake@output[[1]], sep="\t")


