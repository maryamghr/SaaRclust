clusterMapping <- function(soft.clust.probs, cluster.mapping){
	# reading the data table
	dt <- fread(paste("zcat", soft.clust.probs))

	dt.per.clust <- dt[PBchrom %in% paste0("chr", c(1:22, "X")) & PBflag %in% c(0,16),  # keep only PB reads with defined chromosome names and flag 0 or 16
           lapply(.SD, mean),                                               		    # compute the average of (inferred clusters) probabilities per column in each group
           by = .(PBchrom, PBflag),                                         		    # group the data table by original clusters (columns PBchrom and PB flag)
           .SDcols=paste0("V", 1:47)]                                       		    # do these computations only in the inferred clusters (columns names V1 to V47)

	# define the original cluster names
	original.cluster <- paste0(dt.per.clust$PBchrom, "_", dt.per.clust$PBflag)
	# most probable (inferred clusters) for each original cluster
	inferred.cluster <- paste0("V", apply(dt.per.clust[, 3:49], 1, which.max))

	# save the original to inferred chorm mapping in a data table
	write.table(data.table(original.cluster, inferred.cluster), file=cluster.mapping, sep="\t", quote=F)
}


