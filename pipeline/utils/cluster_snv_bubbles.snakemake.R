log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')


library(data.table)

bubbles.cov <- lapply(snakemake@input[["bubble_ss_count"]], fread)

# getting cluster names from the file names
clusters <- sapply(snakemake@input[["bubble_ss_count"]], function(s) strsplit(s, "\\_cluster|\\_snv")[[1]][2])

bubble.clust.cov <- data.table()

for (i in 1:length(bubbles.cov))
{
	colnames(bubbles.cov[[i]]) = c(clusters[i], "bubble")
	if (i==1){
		d <- bubbles.cov[[i]]
	} else {
		d <- merge(d, bubbles.cov[[i]], by="bubble", all=T)
	}
}

# replace NA values by 0
d[is.na(d)] <- 0

d[, max.cov:=max(as.numeric(.SD)), by=bubble]
d[, clust:=clusters[which.max(as.numeric(.SD))], by=.(bubble, max.cov)]
d[, num.max.clust:=length(which(as.numeric(.SD)==max(as.numeric(.SD)))), by=.(bubble, max.cov, clust)]


# calling the cluster with max coverage
d[, bubble.id:=strsplit(bubble, "_")[[1]][1], by=bubble]

# FIXME: properly fix this problem later: there are some multiply mapped bubble chains and some multiallelic bubbles
d[, num.bubble.rep:=nrow(.SD), by=bubble.id]
d <- d[num.bubble.rep<2]
d[, num.bubble.rep:=NULL]

# TODO: replace the above code block by the following line (the above mentioned comment doesn't seem to be a problem at all)
d <- d[, sapply(sum, .SD), by=bubble.id]


# TODO: computing the accuracy of clustering bubbles...
clust.to.chrom <- fread(snakemake@input[["clust_to_chrom"]], header=T)
colnames(clust.to.chrom)[2] <- "clust"
clust.to.chrom[, original.cluster:=strsplit(original.cluster, "_")[[1]][1], by=clust]

d <- merge(d, clust.to.chrom, by="clust")
d[, orig.chrom:=strsplit(bubble, "_")[[1]][2], by=1:nrow(d)]

print("clustering accuracy in the set of bubbles with only one ML cluster:")
print(d[num.max.clust==1, .N, by=original.cluster==orig.chrom])

fwrite(d[num.max.clust==1], file=snakemake@output[["bubbles_clust"]], sep="\t")
