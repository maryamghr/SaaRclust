log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')


library(data.table)

bubbles.cov <- lapply(snakemake@input[["bubble_clust_count"]], function(x) fread(x, col.names=c("N", "clust.forward", "bubble.num", "allele.num", "bubble.chrom", "bubble.flag", "is.reverse")))
clust.partners <- fread(snakemake@input[["clust_to_chrom"]])

# row bind all input bubble cov files (from different libraries and sum up the coverages for every bubble (for both alleles))
d <- Reduce(rbind, bubbles.cov)[, lapply(.SD, sum), .(clust.forward, bubble.num, bubble.chrom, bubble.flag, is.reverse)]

# filter out the bubbles with more than one chrom_flag 
d[, num.rep:=length(unique(paste0(bubble.chrom, "_", bubble.flag))), bubble.num]
d <- d[num.rep == 1][, num.rep:=NULL]

# merge d with clust.partners
d <- merge(d, clust.partners, by="clust.forward")

# switch the cluster of bubbles to the paired cluster if the mapping direction (to SS read) is reverse
d[is.reverse==TRUE, clust.forward:=clust.backward]
d[, `:=`(original.chrom=NULL, clust.backward=NULL)]
d <- merge(d, clust.partners, by="clust.forward")

# sum up coverages for the same bubble/clusters
d[, N:=sum(N), .(clust.forward, bubble.num)] #, bubble.chrom, bubble.flag, original.chrom, clust.backward
d <- d[, head(.SD, 1), .(clust.forward, bubble.num)]

# calling cluster with max coverage
d[, max.clust.cov:=max(N), bubble.num]
d <- d[N==max.clust.cov]

# remove bubbles with more than one max coverage clusters
d[, num.max.cov.clust:=.N, bubble.num]
d <- d[num.max.cov.clust == 1]

# evaluation
print("accuracy among the clustered bubbles (with only one cluster with maximum coverage)")
print(d[, .N, paste0(bubble.chrom, "_", bubble.flag)==original.chrom])

print("accuracy among the clustered bubbles (with only one cluster with maximum coverage per chromosome)")
print(d[, .N, .(paste0(bubble.chrom, "_", bubble.flag)==original.chrom, original.chrom)])

fwrite(d[, .(clust.forward, bubble.num)], file=snakemake@output[["bubbles_clust"]], sep="\t")
