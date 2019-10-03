log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(reshape2)

print('snakemake@input[["clust_to_chrom"]]:')
print(snakemake@input[["clust_to_chrom"]])
print('snakemake@input[["valid_maps"]]:')
print(snakemake@input[["valid_maps"]])
print('snakemake@input[["bubbles_clust"]]')
print(snakemake@input[["bubbles_clust"]])
print('snakemake@output:')
print(snakemake@output)

clust.partners <- fread(snakemake@input[["clust_to_chrom"]])
all.maps <- lapply(snakemake@input[["valid_maps"]], fread)
map <- Reduce(rbind, all.maps)

bubbles.clust <- fread(snakemake@input[["bubbles_clust"]])
colnames(bubbles.clust) <- c("bubbleClust", "bubbleName")

# merge map and bubbles clusters
map <- merge(map, bubbles.clust, by="bubbleName")

# compute SS cluster partners
clust.partners[, SSclust:=clust.forward]
map <- merge(map, clust.partners[, .(SSclust, clust.backward)], by="SSclust")

# keep only the rows in which SS and bubble chroms are the same
map <- map[clust.backward==bubbleClust | SSclust==bubbleClust]

# for each bubble, check whether both alleles are covered by strand seq reads from the same cluster
map[, num.clust.bubble.diff.alleles:=length(unique(bubbleAllele)), .(SSlib, SSclust, bubbleName)]

# have only one row per bubble/SSclust
map <- map[, head(.SD, 1), .(SSlib, SSclust, bubbleName)]

# set bubbleAllele equal to 2 if both alleles of the bubble are covered by the clust
map[num.clust.bubble.diff.alleles==1, num.clust.bubble.diff.alleles:=0]
map[, bubbleAllele:=max(bubbleAllele, num.clust.bubble.diff.alleles), .(SSlib, SSclust, bubbleName)]

# keep a subset of columns
map <- map[, .(SSlib, SSclust, bubbleName, bubbleAllele, clust.backward)]

# convert long to wide data table (put different ss libs in columns)
map <- data.table::dcast(map, SSclust+bubbleName+clust.backward~SSlib, value.var="bubbleAllele")

map[is.na(map)] <- "-"

# split the data table by cluster pairs (chromosome)
map[, first.clust.pair:=min(SSclust, clust.backward), by=SSclust]
map.sp <- split(map, map$first.clust.pair)

outputs <- unlist(snakemake@output)
print(snakemake@output)
print(outputs)


for (d in map.sp){
	print('unique(d$SSclust):')
	print(unique(d$SSclust))
	# split d by SSclust
	d.sp <- split(d, d$SSclust)
	
	lapply(d.sp, function(x) x[, `:=`(SSclust=NULL, clust.backward=NULL, first.clust.pair=NULL)])

	# FIXME: to be rempved later... The problem of not having properly paired clusters should be fixed later on
	if (length(d.sp) != 2){
		print(paste("warning: the size of the cluster pair is", length(d.sp)))
		next()
	}
	
	# make both clusters have the same set of bubbleNames
	d.sp[[1]] <- merge(d.sp[[1]], d.sp[[2]][, .SD, .SDcols="bubbleName"], by="bubbleName", all=T)
	d.sp[[2]] <- merge(d.sp[[2]], d.sp[[1]][, .SD, .SDcols="bubbleName"], by="bubbleName", all=T)

	d.sp[[1]][is.na(d.sp[[1]])] <- "-"
	d.sp[[2]][is.na(d.sp[[2]])] <- "-"

	lapply(d.sp, function(x) setkey(x, bubbleName))
	
	# find the right output file name
	for (i in 1:2)
	{
		cl <- grep(paste0("_cluster", names(d.sp)[i], "_"), outputs)
		fwrite(d.sp[[i]], file=outputs[cl], sep="\t")
	}
}

for (f in outputs)
	system(paste("touch", f))
