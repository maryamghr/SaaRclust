log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(R.utils)

soft.clust.file <- snakemake@input[["soft_clust_file"]]
minimap.file <- snakemake@input[["minimap_file"]]
output.map.file.gz <- snakemake@output[[1]]
output.map.file <- gsub('.gz', '', output.map.file.gz)

load(soft.clust.file)
map <- fread(paste("zcat", minimap.file))
soft.probs <- data.table(soft.clust.obj$soft.pVal)
soft.probs[, PBreadNames := rownames(soft.clust.obj$soft.pVal)]
soft.probs[, ML_cluster_idx := which.max(.SD), by=PBreadNames]
soft.probs[, max_cluster_prob := as.numeric(.SD)[ML_cluster_idx], by=PBreadNames]
soft.probs[, ML_cluster_idx := paste0("V", ML_cluster_idx)]
soft.probs <- soft.probs[, .(PBreadNames, ML_cluster_idx, max_cluster_prob)]

map <- merge(map, soft.probs, by="PBreadNames", allow.cartesian=T)

fwrite(map, file = output.map.file, quote = F, sep = "\t")
gzip(output.map.file, destname=output.map.file.gz)
