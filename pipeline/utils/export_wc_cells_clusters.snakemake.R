log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')


library(data.table)

soft.clust <- get(load(snakemake@input[["soft_clust"]]))
clust.pairs <- fread(snakemake@input[["clust_patrners"]])
min.theta.wc <- snakemake@params[["min_theta_wc"]]

d <- data.table()

for (j in 1:length(soft.clust$theta.param))
{
	cell.theta.wc <- soft.clust$theta.param[[j]][,3]
	wc.clusters <- which(cell.theta.wc > min.theta.wc)
	
	if (length(wc.clusters) == 0) next

	d <- rbind(d, data.table(lib=names(soft.clust$theta.param)[j], thetawc=cell.theta.wc[wc.clusters], clust.forward=paste0('V', wc.clusters)))
}

clust.pairs[, first.pair:=min(clust.forward, clust.backward), by=clust.forward]

d <- merge(d, clust.pairs, by='clust.forward')

# compute the number of wc clusters per clust pair for every lib
d[, num_wc_clust_per_clust_pair:=.N, by=.(lib, first.pair)]

# take only the cell/clusters in which both cluster pairs are wc, and subset the required columns
d <- d[num_wc_clust_per_clust_pair==2, .(lib, thetawc, clust.forward, clust.backward)]

fwrite(d, file=snakemake@output[[1]], sep='\t')
