args=commandArgs(TRUE)

.libPaths( c( .libPaths(), args[4]))

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

### finding the garbage cluster
soft.clust <- get(load(args[2]))
# getting cluster partners
partners <- findClusterPartners(soft.clust$theta.param)
# finding the unpaored cluster
garbage.clust <- paste0("V", setdiff(1:(length(partners)+1), partners))

counts <- lapply(args[1], fread)

valid.chr_flag.names <- paste0(paste0("chr", c(1:22, "X")), "_", c(rep(0,23), rep(16, 23)))

for (i in 1:length(counts))
{
	colnames(counts[[i]]) <- c("count", "original.chrom", "clust.forward")
	counts[[i]] <- counts[[i]][original.chrom %in% valid.chr_flag.names & clust.forward != garbage.clust]
}

total.counts <- Reduce(rbind, counts)[, lapply(.SD, sum), by=.(original.chrom, clust.forward)]

clust.to.chr_dir.map <- total.counts[, rank:=rank(-count), by=clust.forward][rank==1]

# adding chromosome to the columns
clust.to.chr_dir.map[, chrom:=strsplit(original.chrom, "_")[[1]][1], by=original.chrom]

# finding cluster pairs
clust.to.chr_dir.map[, clust.backward:=rev(clust.forward), by=chrom]

fwrite(clust.to.chr_dir.map[, .(original.chrom, clust.forward, clust.backward)], args[3], sep="\t")
