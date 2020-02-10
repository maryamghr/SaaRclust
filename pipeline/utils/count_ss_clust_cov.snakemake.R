log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source('utils/count_short_read_clust_cov.R')

minimap <- fread(paste('zcat', snakemake@input[['minimap']]), select=c('SSreadNames', 'strand', 'PBreadNames'))
soft.clust <- get(load(snakemake@input[['soft_clust']]))
clust.partners <- fread(snakemake@input[['clust_partners']])

ss.clust.cov <- getClustCovInShortReads(minimap, soft.clust, clust.partners)

fwrite(ss.clust.cov, file=snakemake@output[[1]], sep='\t', row.names=F)

