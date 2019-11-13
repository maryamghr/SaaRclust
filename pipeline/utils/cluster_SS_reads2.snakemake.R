log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source("utils/clustering_functions.R")

print('reading the ss clusters count')
cov <- lapply(snakemake@input[["ss_clust_count"]], fread)
print('reading the cluster partners')
clust.partners <- fread(snakemake@input[["clust_partners"]])

print('merging coverages')
d <- mergeAllCoverages(cov, clust.partners)
print('calling majority vote clusters')
d <- callMajorityVoteCluster(d)
print('evaluating clustering')
evaluateSSclustering(d)

print('writing the output')
fwrite(d[, .(clust.forward, name)], file=snakemake@output[[1]], sep="\t")


