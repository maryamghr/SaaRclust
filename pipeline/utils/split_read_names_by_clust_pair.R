log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(dplyr)

soft.clust.filename <- snakemake@input[['soft_clust']]
#dir.name <- dirname(soft.clust.filename)
dir.name <- dirname(snakemake@output[[1]])

clust.pairs <- fread(snakemake@input[['clust_pairs']])
clust.pairs[, clust_pair:=paste0(min(clust.forward, clust.backward), '_', max(clust.forward, clust.backward)), by=1:nrow(clust.pairs)]
soft.clust <- fread(soft.clust.filename)
#soft.clust[, haplo:='H1']

if (snakemake@params[["type"]] == "long_read"){
  colnames(soft.clust) = c('read_name', 'clust.forward', 'prob')
}
if (snakemake@params[["type"]] %in% c("bubble", "SS")){
  colnames(soft.clust)[1:2] = c('clust.forward', 'read_name')
}

soft.clust <- merge(soft.clust, clust.pairs[, .(clust.forward, clust_pair)], by='clust.forward')
soft.clust.sp <- split(soft.clust, soft.clust[, clust_pair])

#lapply(1:length(soft.clust.sp), function(x) fwrite(soft.clust.sp[[x]][, .(read_name, haplo,clust.forward,prob)], file=file.path(dir.name, paste0('cluster', names(soft.clust.sp)[x], '.data')), sep="\t", col.names=F)) %>% invisible
lapply(1:length(soft.clust.sp), function(x) fwrite(soft.clust.sp[[x]][, .(read_name,clust.forward)], file=file.path(dir.name, paste0('cluster', names(soft.clust.sp)[x], '.data')), sep="\t", col.names=F)) %>% invisible
