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
print('snakemake@input[["ss_clust"]]')
print(snakemake@input[["ss_clust"]])
print('snakemake@output:')
print(snakemake@output)

wc.cell.clust <- fread(snakemake@input[["wc_cell_clust"]])
clust.partners <- fread(snakemake@input[["clust_to_chrom"]])
all.maps <- lapply(snakemake@input[["valid_maps"]], fread)
map <- Reduce(rbind, all.maps)

bubbles.clust <- fread(snakemake@input[["bubbles_clust"]])
# add bubble id...
bubbles.clust[, bubble.id:=sapply(name, function(x) as.integer(strsplit(x, '_')[[1]][2]))]
bubbles.clust <- bubbles.clust[, head(.SD, 1), by=bubble.id]
bubbles.clust <- bubbles.clust[, .(clust.forward, bubble.id)]
colnames(bubbles.clust) <- c("bubbleClust", "bubbleName")

ss.clust <- fread(snakemake@input[["ss_clust"]])
colnames(ss.clust) <- c("SSclust", "SSname")

# merge map and bubbles clusters
map <- merge(map, bubbles.clust, by="bubbleName")
map <- merge(map, ss.clust, by="SSname")

# compute SS cluster partners
clust.partners[, SSclust:=clust.forward]
map <- merge(map, clust.partners[, .(SSclust, clust.backward)], by="SSclust")

# keep only the rows in which SS and bubble chroms are the same
map <- map[clust.backward==bubbleClust | SSclust==bubbleClust]

# for each bubble, check whether both alleles are covered by strand seq reads from the same cluster
map[, num.clust.bubble.diff.alleles:=length(unique(bubbleAllele)), .(SSlib, SSclust, bubbleName)]

# have only one row per bubble/SSclust/SSlib
map <- map[, head(.SD, 1), .(SSlib, SSclust, bubbleName)]

# set bubbleAllele equal to 2 if both alleles of the bubble are covered by the clust
map[num.clust.bubble.diff.alleles>1, bubbleAllele:=2]

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
	clust.pair = unique(d$SSclust)
	print(clust.pair)

	# split d by SSclust
	d.sp <- split(d, d$SSclust)
	
	lapply(d.sp, function(x) x[, `:=`(SSclust=NULL, clust.backward=NULL, first.clust.pair=NULL)])

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

	# select a subset of cells that are wc in this cluster pair
	wc.cells <- wc.cell.clust[clust.forward==clust.pair[1], lib]
	selected.col.names <- c('bubbleName', wc.cells)

	d.sp[[1]] <- d.sp[[1]][, selected.col.names, with=FALSE]
	d.sp[[2]] <- d.sp[[2]][, selected.col.names, with=FALSE]


	
	# find the right output file name
	for (i in 1:2)
	{
		cl <- grep(paste0("cluster", names(d.sp)[i], "_"), outputs)
		print(paste('cl =', cl, ', file =', outputs[cl]))
		fwrite(d.sp[[i]], file=outputs[cl], sep="\t")
	}
}

# touching the garbage cluster
garbage.clust <- setdiff(paste0('V',1:(nrow(clust.partners)+1)), clust.partners[, clust.forward])
garbage.out.file.idx <- grep(paste0("cluster", garbage.clust, "_"), outputs)
system(paste("touch", outputs[garbage.out.file.idx]))
