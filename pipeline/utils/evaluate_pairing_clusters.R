args=commandArgs(TRUE)

.libPaths( c( .libPaths(), args[4]))

suppressPackageStartupMessages(library(SaaRclust))
suppressPackageStartupMessages(library(lpSolve))
suppressPackageStartupMessages(library(data.table))


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


numFoundClusters <- function (ord, chr, flag) 
{
    hard.clust.dt <- data.table(name = names(ord), clust = ord, chrom = chr, flag = flag)
    hard.clust.dt[, `:=`(chrom_flag_count, .N), by = .(clust, chrom, flag)]
    hard.clust.to.chrom <- hard.clust.dt[, head(.SD, 1), by = .(clust, chrom, flag)]
    hard.clust.to.chrom[, name:=NULL]
    hard.clust.to.chrom[, `:=`(chrom_flag_rank, rank(-chrom_flag_count)), by = clust]
    hard.clust.to.chrom <- hard.clust.to.chrom[chrom_flag_rank == 1]
    hard.clust.to.chrom.unq <- unique(hard.clust.to.chrom[, .(chrom, flag)])
    return(nrow(hard.clust.to.chrom.unq))
}


findClusterPartners <- function (theta.param = NULL) 
{
    num.clusters <- nrow(theta.param[[1]])
    if (num.clusters%%2 != 0) {
        theta.sums <- Reduce("+", theta.param)
        remove.clust <- which.max(theta.sums[, 3])
        message("    Removed cluster ", remove.clust, " to ensure even number of clusters!!!")
        theta.param <- lapply(theta.param, function(x) x[-remove.clust, 
            ])
    }
    theta.param.wc <- lapply(theta.param, function(x) x[, 3])
    all.theta.param.wc <- do.call(cbind, theta.param.wc)
    d <- as.matrix(dist(all.theta.param.wc))
    d <- max(d) - d
    diag(d) <- 0
    max.partners <- lpSolve::lp.assign(d, "max")
    max.partners.m <- max.partners$solution
    max.partners.idx <- which(max.partners.m > 0, arr.ind = TRUE)
    max.partners.idx <- max.partners.idx[max.partners.idx[, 1] < 
        max.partners.idx[, 2], ]
    colnames(max.partners.idx) <- c("Cluster1", "Cluster2")
    return(max.partners.idx)
}

