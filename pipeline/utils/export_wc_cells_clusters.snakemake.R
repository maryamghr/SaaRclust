log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)

soft.clust <- get(load(snakemake@input[[1]]))

d <- data.table()

for (j in 1:length(soft.clust$theta.param))
{
	wc.clusters <- which(exp(soft.clust$theta.param[[j]][,3])>0.9)
	wc.clusters <- paste0('V', wc.clusters)
	d <- rbind(d, data.table(cell=names(soft.clust$theta.param)[j], cluster=wc.clusters))
}

fwrite(d, file=snakemake@output[[1]], sep='\t')
