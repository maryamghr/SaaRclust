library('data.table')

getClustCovInShortReads <- function(minimap, soft.clust, clust.partners, out.file){
	start.time = Sys.time()
	print('computing the ML clusters')
	long.read.names <- rownames(soft.clust$soft.pVal)
	soft.clust <- data.table(soft.clust$soft.pVal)
	soft.clust[, clust.forward := which.max(.SD), by=1:nrow(soft.clust)]
	soft.clust[, clust.prob := as.numeric(.SD)[clust.forward], by=1:nrow(soft.clust)]
	soft.clust[, clust.forward := paste0("V", clust.forward)]
	soft.clust[, long.read.names:=long.read.names]
	soft.clust <- soft.clust[, .(long.read.names, clust.forward, clust.prob)]

	print(Sys.time()-start.time)
	start.time = Sys.time()
	
	print('counting the clusters coverage in short reads')

	colnames(minimap) <- c("short.read.names", "strand", "long.read.names")
	
	short.read.clust <- merge(minimap, soft.clust, by='long.read.names')
	short.read.clust <- merge(short.read.clust, clust.partners, by='clust.forward')
	short.read.clust[strand=='-', clust.forward:=clust.backward]

	# computing the cluster counts in each short read
	short.read.clust[, clust.cov:=.N, by=.(short.read.names, clust.forward)]

	# removing repetitive read/clust rows in the data table
	short.read.clust <- short.read.clust[, head(.SD, 1), by=.(short.read.names, clust.forward)]
	print(Sys.time()-start.time)

	return(short.read.clust[, .(short.read.names, clust.forward, clust.cov)])
}
