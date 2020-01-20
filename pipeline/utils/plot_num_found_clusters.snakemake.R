log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), 'utils/R-packages/'))

suppressPackageStartupMessages(library(SaaRclust))
library(data.table)

hard.clusters <- lapply(snakemake@input, function(x) get(load(x)))

dt <- data.table()

for (i in 1:length(hard.clusters)){
	filename.sp <- strsplit(snakemake@input[[i]], '_')
	num.clust <- as.numeric(strsplit(filename.sp[[1]][length(filename.sp[[1]])], 'clusters.RData')[[1]][1])
	num.clust.unmerged.found <- numFoundClusters(hard.clusters[[i]]$unmerged.ord, hard.clusters[[i]]$pb.chr, hard.clusters[[i]]$pb.flag)
	num.clust.found <- numFoundClusters(hard.clusters[[i]]$ord, hard.clusters[[i]]$pb.chr, hard.clusters[[i]]$pb.flag)
	dt <- rbind.data.frame(dt, data.table(num.clust=num.clust, num.clust.unmerged.found=num.clust.unmerged.found, num.clust.found=num.clust.found))
}

dt.long <- melt(dt, id.var='num.clust', variable.name='merged', value.name='num.found.clust')
dt.long[merged=='num.clust.unmerged.found', merged:='not merged']
dt.long[merged=='num.clust.found', merged:='merged']

plt <- ggplot(dt.long, aes(x=num.clust, y=num.found.clust, color=merged)) + geom_point() + geom_line()

# creating output files
fwrite(dt, file=snakemake@output[["table"]], sep="\t")
ggsave(filename=snakemake@output[["plot"]], plt)
