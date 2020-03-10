log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(R.utils)

soft.clust.file <- snakemake@input[["soft_clust_file"]]

soft.clust <- get(load(soft.clust.file))

# keep.rowname = TRUE will retain the rownames of that object in a column named rn. So, here the "rn" column shows the long read names.
start.time = Sys.time()
print('computing the ML clusters')
long.read.names <- rownames(soft.clust$soft.pVal)
soft.clust <- data.table(soft.clust$soft.pVal)
soft.clust[, clust.forward := which.max(.SD), by=1:nrow(soft.clust)]
soft.clust[, clust.prob := as.numeric(.SD)[clust.forward], by=1:nrow(soft.clust)]
soft.clust[, clust.forward := paste0("V", clust.forward)]
soft.clust[, long.read.names:=long.read.names]
soft.clust <- soft.clust[, .(long.read.names, clust.forward, clust.prob)]


fwrite(soft.clust, file=snakemake@output[[1]], sep="\t", quote = F, row.names = F, col.names = F)
