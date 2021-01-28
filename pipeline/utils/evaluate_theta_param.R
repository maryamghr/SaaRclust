library(data.table)
library(ggplot2)
library(gridExtra)

countfiles.dir <- '../../HG00733/SaaRclust/ground_truth_strand_states/'
soft.clust <- get(load('/MMCI/TM/scratch/maryam/haploclust-css-HG00514/SaaRclust/pipeline/aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00514/Clusters/HG00514_chunk000_clusters.RData'))
hard.clust <- get(load('/MMCI/TM/scratch/maryam/haploclust-css-HG00514/SaaRclust/pipeline/aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00514/Clusters/hardClusteringResults_100clusters.RData'))
clust.to.chrom <- fread('../../HG00733/SaaRclust/Clusters/clust_partners.txt')
clust.to.chrom[,cluster:=clust.forward]

soft.theta <- soft.clust$theta.param
hard.theta <- hard.clust$theta.param

compare_theta_with_wfrac <- function(theta, clust.to.chrom, countfiles.dir="", title="", num.clusters=80)
{
	theta.wc = data.table()

	for (cell in 1:length(theta))
	{
		libname = names(theta)[cell]
	#	cell.theta.wc = data.table(lib=libname, clust=paste0('V', 1:47), thetawc=theta[[cell]][,3])
		cell.theta.wc = data.table(lib=libname, clust=1:num.clusters, thetawc=theta[[cell]][,3])
		cell.theta.wc = merge(cell.theta.wc, clust.to.chrom[, .(chrom, clust)], by="clust")
	#	cell.theta.wc[, chrom:=sapply(original.chrom, function(x) strsplit(x, '_')[[1]][1])]
		
		# reading W and C read counts data
		countfile = paste0(countfiles.dir, libname, '_chrom_haplo_count.data')
		counts = fread(countfile)

		cell.theta.wc <- merge(cell.theta.wc, counts[, .(chrom, W.frac=W/(W+C))], by="chrom")

		theta.wc = rbind(theta.wc, cell.theta.wc[, .(lib, chrom, clust, thetawc, W.frac)])
	}

	return(list(theta.wc, ggplot(theta.wc, aes(x=W.frac, y=thetawc))+geom_point()+ggtitle(title)))
}

soft.clust.plt <- compare_theta_with_wfrac(soft.theta, 'soft clustering theta')[[2]]
hard.clust.plt <- compare_theta_with_wfrac(hard.theta, 'hard clustering theta')[[2]]

grid.arrange(soft.clust.plt, hard.clust.plt, nrow=2)


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


# plotting number of found clusters
plot_num_found_clusters <- function(hard.clusters, num.clusters)
{
	dt <- data.table()

	for (i in 1:length(hard.clusters)){
		num.clust <- num.clusters[i]
		num.clust.unmerged.found <- numFoundClusters(hard.clusters[[i]]$unmerged.ord, hard.clusters[[i]]$pb.chr, hard.clusters[[i]]$pb.flag)[[2]]
		num.clust.found <- numFoundClusters(hard.clusters[[i]]$ord, hard.clusters[[i]]$pb.chr, hard.clusters[[i]]$pb.flag)[[2]]
		dt <- rbind.data.frame(dt, data.table(num.clust=num.clust, num.clust.unmerged.found=num.clust.unmerged.found, num.clust.found=num.clust.found))
	}

	dt.long <- melt(dt, id.var='num.clust', variable.name='merged', value.name='num.found.clust')
	dt.long[merged=='num.clust.unmerged.found', merged:='not merged']
	dt.long[merged=='num.clust.found', merged:='merged']

	plt <- ggplot(dt.long, aes(x=num.clust, y=num.found.clust, color=merged)) + geom_point() + geom_line()
	
	return(plt)
}

findClusterPartners <- function(theta.param = NULL) 
{
	num.clusters <- nrow(theta.param[[1]])
	rownames(theta.param[[1]]) <- 1:num.clusters
	if (num.clusters%%2 != 0) {
	theta.sums <- Reduce("+", theta.param)
	remove.clust <- which.max(theta.sums[, 3])
	message("    Removed cluster ", remove.clust, " to ensure even number of clusters!!!")
	theta.param <- lapply(theta.param, function(x) x[-remove.clust,])
	}
	theta.param.wc <- lapply(theta.param, function(x) x[, 3])
	all.theta.param.wc <- do.call(cbind, theta.param.wc)
	d <- as.matrix(dist(all.theta.param.wc))
	d <- max(d) - d
	diag(d) <- 0
	max.partners <- lpSolve::lp.assign(d, "max")
	max.partners.m <- max.partners$solution
	max.partners.idx <- which(max.partners.m > 0, arr.ind = TRUE)
	max.partners.idx <- max.partners.idx[max.partners.idx[, 1] < max.partners.idx[, 2], ]

	get.clust.idx <- as.numeric(rownames(theta.param[[1]]))
	colnames(max.partners.idx) <- c("Cluster1", "Cluster2")
	return(matrix(get.clust.idx[max.partners.idx], ncol=2))
}


num.clusters.plt <- plot_num_found_clusters(list(hard.clust), 100)

partners <- findClusterPartners(soft.clust$theta.param)

# making a data table of cluster partners
clust.partners <- data.table(clust=c(partners[,1], partners[,2]), pair=c(partners[,2], partners[,1]))

hard.clust.to.chromflag <- numFoundClusters(hard.clust$ord, hard.clust$pb.chr, hard.clust$pb.flag)[[1]]
hard.clust.to.chromflag <- merge(hard.clust.to.chromflag, clust.partners, by="clust", all=TRUE)

hard.clust.to.chromflag[, `:=`(original.chrom=paste0(chrom, '_', flag), clust.forward=paste0('V',clust), clust.backward=paste0('V',pair))]





