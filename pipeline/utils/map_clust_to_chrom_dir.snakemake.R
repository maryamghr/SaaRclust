args=commandArgs(TRUE)

.libPaths( c( .libPaths(), args[3]))

suppressPackageStartupMessages(library(SaaRclust))
suppressPackageStartupMessages(library(lpSolve))
suppressPackageStartupMessages(library(data.table))


## new version of the function after fixing the bug

findClusterPartners <- function(theta.param=NULL) {
    
	## If there is an uneven number of clusters remove the one with the most WC states
	num.clusters <- nrow(theta.param[[1]])
	rownames(theta.param[[1]]) <- 1:num.clusters

	if (num.clusters %% 2 != 0) {
		#Find cluster with WC state in majority of cells
		theta.sums <- Reduce("+", theta.param)
		remove.clust <- which.max(theta.sums[,3])
		message("    Removed cluster ", remove.clust, " to ensure even number of clusters!!!")
		theta.param <- lapply(theta.param, function(x) x[-remove.clust,])
	}

	## get only wc thetas
	theta.param.wc <- lapply(theta.param, function(x) x[,3])
	## cbind wc thetas for all single cells
	all.theta.param.wc <- do.call(cbind, theta.param.wc)
	## compute the pairwise distance of all clusters wc thetas
	d <- as.matrix(dist(all.theta.param.wc))
	## convert distance to a similarity measure
	d <- max(d) - d
	## set diagonal values to zero
	diag(d) <- 0
	## Find pairs of clusters with the highest similarity
	max.partners <- lpSolve::lp.assign(d, "max")
	max.partners.m <- max.partners$solution #matrix with pairs of clusters with maximal similarity
	## Extract indices of pair of clusters
	max.partners.idx <- which(max.partners.m > 0, arr.ind = TRUE)
	max.partners.idx <- max.partners.idx[max.partners.idx[,1] < max.partners.idx[,2],] #remove duplicate cluster partners
	
	get.clust.idx <- as.numeric(rownames(theta.param[[1]]))
	colnames(max.partners.idx) <- c("Cluster1", "Cluster2")
	return(matrix(get.clust.idx[max.partners.idx], ncol=2))
}

numFoundClusters <- function (ord, chr, flag) 
{
	hard.clust.dt <- data.table(name = names(ord), clust = ord, chrom = chr, flag = flag)
	hard.clust.dt[, `:=`(chrom_flag_count, .N), by = .(clust, chrom, flag)]
	hard.clust.to.chrom <- hard.clust.dt[, head(.SD, 1), by = .(clust, chrom, flag)]
	hard.clust.to.chrom[, name:=NULL]
	hard.clust.to.chrom[, `:=`(chrom_flag_rank, rank(-chrom_flag_count)), by = clust]
	hard.clust.to.chrom <- hard.clust.to.chrom[chrom_flag_rank == 1]
	hard.clust.to.chromflag <- hard.clust.to.chrom[,.(clust, chrom, flag)][order(chrom)]
	hard.clust.to.chrom.unq <- unique(hard.clust.to.chrom[, .(chrom, flag)])
	return(list(hard.clust.to.chromflag, nrow(hard.clust.to.chrom.unq)))
}


hard.clust <- get(load(args[1]))
# getting cluster partners
partners <- findClusterPartners(hard.clust$theta.param)

clust.partners <- data.table(clust=c(partners[,1], partners[,2]), pair=c(partners[,2], partners[,1]))

hard.clust.to.chromflag <- numFoundClusters(hard.clust$ord, hard.clust$pb.chr, hard.clust$pb.flag)[[1]]
hard.clust.to.chromflag <- merge(hard.clust.to.chromflag, clust.partners, by="clust", all=TRUE)

hard.clust.to.chromflag[, `:=`(original.chrom=paste0(chrom, '_', flag), clust.forward=paste0('V',clust), clust.backward=paste0('V',pair))]

fwrite(hard.clust.to.chromflag[, .(original.chrom, clust.forward, clust.backward)], args[2], sep="\t")
