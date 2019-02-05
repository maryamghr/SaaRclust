log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source("utils/clustering_functions.R")

cov <- lapply(snakemake@input[["ss_clust_count"]], fread)
clust.to.chrom <- fread(snakemake@input[["clust_to_chrom"]])
clust.partners <- fread(snakemake@input[["clust_partners"]])[, 2:3]
clust.partners.rev.col <- clust.partners
colnames(clust.partners.rev.col) <- colnames(clust.partners)[2:1]
clust.partners <- rbind(clust.partners, clust.partners.rev.col)
colnames(clust.partners) <- c("clust", "pair")

d <- mergeAllCoverages(cov, clust.partners)
d <- callMajorityVoteCluster(d)
evaluateSSclustering(d, clust.to.chrom)

fwrite(d[num.max.clust==1], file=snakemake@output[["bubbles_clust"]], sep="\t")
