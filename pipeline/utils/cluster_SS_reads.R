clusterSSreads <- function(map, clust.to.chrom, clust.pairs) {
	# create a dataframe listing each cluster in the first column with its paired cluster in the second column
	cluster.pairs <- data.table(pair1=c(clust.pairs$clust.forward, clust.pairs$clust.backward),
		                    pair2=c(clust.pairs$clust.backward, clust.pairs$clust.forward))

	# find the garbage (unpaired cluster)
	garbage.clust <- which(is.na(match(paste0("V", 1:47), cluster.pairs$pair1)))
	garbage.clust <- paste0("V", garbage.clust)


	# wierd behaviour in reading the first line!!!
	# Warning message:
	#In fread("/local/home/maryam/research/saarclust/haploclust/aln_NW130711_291_sorted.maf") :
	#  Detected 1 column names but the data has 16 columns (i.e. invalid file). Added 15 extra default column names at the end.

	names(map) <- strsplit(names(map)[1], split = "\t")[[1]]

	# sort map by SSreadNames (check later in the pipeline why it is not sorted)
	setkey(map, SSreadNames)

	# kick out SS reads with unknown chromosomes
	map <- map[SSchrom %in% paste0("chr", c(1:22, "X"))]

	# convert SSflag column to 0 and 16
	map[, SSflag:=bitwAnd(SSflag, 16)]

	# add a new column in map for the true PB cluster
	map[, truePBclust:= rep(clust.to.chrom[original.cluster==paste0(c(PBchrom[1], PBflag[1]), collapse="_"), inferred.cluster], .N),
	    by=.(PBchrom, PBflag)]
	# Note that for the PB reads with unknown chromosome or unknown flags (neither 0 nor 16), truePBclust is NA

	# add a new column in map for the true SS cluster
	map[, trueSSclust:= rep(clust.to.chrom[original.cluster==paste0(c(SSchrom[1], SSflag[1]), collapse="_"), inferred.cluster], .N),
	    by=.(SSchrom, SSflag)]

	# add a new column for SS clust (based on strand and the ML clust of the PB read mapped to the SS read)
	# Do this operation for each group of ML_cluster_idx and strand
	map[, SSclust := rep(ifelse(strand[1]=="+" & ML_cluster_idx != garbage.clust, ML_cluster_idx[1], cluster.pairs[pair1==ML_cluster_idx[1], pair2]), .N),
	    by = .(ML_cluster_idx, strand)]
	# Note that for the lines where the ML_clust is the garbage cluster, SSclust is NA

	# shrink the table and assign a score to each cluster (sum of soft clust PB probs for that clust) assigned to every SS
	SS.clust <- map[, c(head(.SD, 1), clustWeight = sum(max_cluster_prob)),
		        by = .(SSreadNames, SSclust),
		        .SDcols = c("SSlibNames", "SSflag", "SSchrom", "SSpos", "trueSSclust")]

	# add a column for rank of SS assigned clusters based on SS cluster weights for each SS read
	SS.clust[, clustRank:= order(clustWeight, decreasing = T),
		        by = SSreadNames]

	# add a new column counting the number of assigned clusters for each SS read
	SS.clust[, numAssignedClust:=.N,
		 by = .(SSreadNames)]

	# cluster assignment to SS reads by subsetting the table to the clusters with first rank
	SS.clust.assign <- SS.clust[clustRank==1]

	# compute the fraction of wrongly assigned SS reads
	# It doesn't count SS reads with NA SSclust
	SS.clust.assign[SSclust != trueSSclust, .N] / SS.clust.assign[, .N]

	# compare SS reads with true and wrong chrom assignment
	# add a new column to SS.clust.assign for the correctness of the cluster assignment:
	SS.clust.assign[, trueClustAssigned:=(SSclust == trueSSclust)]
	SS.clust.assign[, clustRank:=NULL]

	return(SS.clust.assign)
}


