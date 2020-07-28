log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
source('utils/count_short_read_clust_cov.R')


minimap <- fread(paste('zcat', snakemake@input[['minimap']]), select=c(1,5,6), fill=T)
# TODO: correct it later after removing the ground true info from the pipeline
# remove the ground true part of the long read names
minimap[, name:=sapply(V6, function(x) strsplit(x, '/ccs')[[1]][1])]
minimap[, name:=paste0(name, '/ccs')]
minimap <- minimap[, c(1,2,4), with=F]

if (endsWith(snakemake@input[['soft_clust']],".data")){
  soft.clust <- fread(snakemake@input[['soft_clust']])
} else {
  soft.clust <- get(load(snakemake@input[['soft_clust']]))
}

clust.partners <- fread(snakemake@input[['clust_partners']])

bubble.clust.cov <- getClustCovInShortReads(minimap, soft.clust, clust.partners)
# get bubble ids from the full bubble names
#bubble.clust.cov <- bubble.clust.cov[, bubble.id:=sapply(short.read.names, function(x) strsplit(x, '_')[[1]][2])]
#bubble.clust.cov <- bubble.clust.cov[, bubble.id:=short.read.names]

fwrite(bubble.clust.cov[, .(short.read.names, clust.forward, clust.cov)], file=snakemake@output[[1]], sep='\t', row.names=F)

