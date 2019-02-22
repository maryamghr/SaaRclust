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


