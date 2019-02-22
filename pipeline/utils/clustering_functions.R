library(data.table)
library(reshape2)
library(stringr)
library(binaryLogic)


mergeAllCoverages <- function(cov, clust.partners) {
	for (i in 1:length(cov))
	{
		colnames(cov[[i]]) <- c("cov", "name", "lib", "flag", "chrom", "pos", "len", "dir", "clust.forward")
		# adding cluster partners
		cov[[i]] <- merge(cov[[i]], clust.partners, by="clust.forward")
		# replace each cluster by its pair if the mappint direction is minus
		cov[[i]][dir=="-", clust.forward:=clust.backward]
		# keep useful columns
		cov[[i]] <- cov[[i]][, .(name, lib, flag, chrom, pos, len, clust.forward, cov, original.chrom)]
		# summing up the coverages in the repeated name/clust columns
		cov[[i]][, cov:=sum(cov), by=.(name, clust.forward)]
		cov[[i]] <- cov[[i]][, head(.SD, 1), by=.(name, clust.forward)]
		cov[[i]][, original.chrom:=NULL]
		# merge the ground true clusters with correctly oriented cov
		cov[[i]] <- merge(cov[[i]], clust.partners, by="clust.forward")
		cov[[i]][, clust.backward:=NULL]
	}

	return(Reduce(rbind, cov)[, lapply(.SD, sum), by=.(name, lib, flag, chrom, pos, len, clust.forward, original.chrom)])
}


callMajorityVoteCluster <- function(d) {
	# getting the max coverage of all clusters for each SS read
	d[, max.cov:=max(cov), by=name]
	# subsetting only the clusters with max coverage for each SS read
	d <- d[cov==max.cov]
	# counting the number of maximum clusters for each SS read
	d <- d[, num.max.clust:=.N, by=name]
	colnames(d) <- c("name", "lib", "flag", "chrom", "pos", "len", "clust.forward", "inferred.chrom.flag", "cov", "max.cov", "num.max.clust")
	return(d[num.max.clust==1])
}


evaluateSSclustering <- function(d) {
	# compute SS mapping direction based on the flag
	# compute binary flag values
	d[, binary.flag:=paste(as.binary(flag[1]), collapse=""), by=flag]
	# get the forth (direction flag) digit from the right
	d[, dir.flag:=rep(substr(binary.flag[1], nchar(binary.flag[1])-4, nchar(binary.flag[1])-4), .N), by=binary.flag]
	# set the dir flag to 0 for the flags that have less than 5 binaray characters
	d[dir.flag=="", dir.flag:=0]
	# convert dir flag 1 to 16
	d[dir.flag==1, dir.flag:=16]

	d[, name:=paste0(name, "_", lib, "_", flag, "_", chrom, "_", pos, "_len:", len)]

	print("clustering accuracy among SS reads with only one ML cluster:")
	print(d[, .N, by=paste0(chrom, "_", dir.flag)==inferred.chrom.flag])
}
