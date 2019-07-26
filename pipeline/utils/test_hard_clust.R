.libPaths( c( .libPaths(), 'utils/R-packages/'))

suppressPackageStartupMessages(library(SaaRclust))
library(data.table)
library(ggplot2)
library(gridExtra)
library(pheatmap)


getclusttochrom <- function (ord, chr, flag) 
{
    hard.clust.dt <- data.table(name = names(ord), 
        clust = ord, chrom = chr, flag = flag)
    hard.clust.dt[, `:=`(chrom_flag_count, .N), by = .(clust, 
        chrom, flag)]
    hard.clust.to.chrom <- hard.clust.dt[, head(.SD, 1), by = .(clust, 
        chrom, flag)][, `:=`(name, NULL)]
    hard.clust.to.chrom[, `:=`(chrom_flag_rank, rank(-chrom_flag_count)), 
        by = clust]
    hard.clust.to.chrom <- hard.clust.to.chrom[chrom_flag_rank == 
        1]
    return(hard.clust.to.chrom)
}




hard.cluster <- get(load('aligns_k15_w-default_f0.0002_z500/SaaRclust_results_HG00514/Clusters/hardClusteringResults_100clusters.RData'))
best.alignments <- get(load('aligns_k15_w-default_f0.0002_z500/SaaRclust_results_HG00514/RawData/representativeAligns.RData'))

best.alignments$PBreadNames <- factor(best.alignments$PBreadNames, levels=unique(best.alignments$PBreadNames))
    
#split data by Strand-seq library
tab.l <- split(best.alignments, best.alignments$SSlibNames)

### Count directional reads ###
counts.l <- countDirectionalReads(tab.l)


cluster.to.chrom.unmerged <- getclusttochrom(hard.cluster$unmerged.ord, hard.cluster$pb.chr, hard.cluster$pb.flag)
cluster.to.chrom <- getclusttochrom(hard.cluster$ord, hard.cluster$pb.chr, hard.cluster$pb.flag)

theta.l <- estimateTheta(counts.l, hard.clust=hard.cluster$unmerged.ord, alpha=0.05)


theta.all <- do.call(cbind, theta.l)
hc <- hclust(dist(theta.all))
hc.clust <- cutree(hc, k=47)

cluster.to.chrom.unmerged[, merged_clust:=hc.clust[clust]]


## collecting statisrics from SS read counts for PB reads
pb.cov <- rowSums(Reduce(cbind, counts.l))



# digging into theta estims for chr 21
theta.chr21.16 <- Reduce(rbind, lapply(theta.l, function(x) x[56,]))
rownames(theta.chr21.16) <- names(theta.l)
theta.chr21.0 <- Reduce(rbind, lapply(theta.l, function(x) x[57,]))
rownames(theta.chr21.0) <- names(theta.l)


theta.chr21.0.plt <- ggplot(data = melt(theta.chr21.0), aes(x=Var2, y=Var1, fill=value, label=value))+geom_tile(color = "black")+scale_fill_gradient2(low = "red", high = "blue")+ggtitle("theta chr 21 forward")+xlab("strand states")+ylab("cells")
theta.chr21.16.plt <- ggplot(data = melt(theta.chr21.16), aes(x=Var2, y=Var1, fill=value, label=value))+geom_tile(color = "black")+scale_fill_gradient2(low = "red", high = "blue")+ggtitle("theta chr 21 backward")+xlab("strand states")+ylab("cells")

grid.arrange(theta.chr21.0.plt, theta.chr21.16.plt, nrow=1, ncol=2)



theta.chr1.16 <- Reduce(rbind, lapply(theta.l, function(x) x[1,]))
rownames(theta.chr1.16) <- names(theta.l)
theta.chr1.0 <- Reduce(rbind, lapply(theta.l, function(x) x[9,]))
rownames(theta.chr1.0) <- names(theta.l)
theta.chr1.0.plt <- ggplot(data = melt(theta.chr1.0), aes(x=Var2, y=Var1, fill=value, label=value))+geom_tile(color = "black")+scale_fill_gradient2(low = "red", high = "blue")+ggtitle("theta chr 21 forward")+xlab("strand states")+ylab("cells")
theta.chr1.16.plt <- ggplot(data = melt(theta.chr1.16), aes(x=Var2, y=Var1, fill=value, label=value))+geom_tile(color = "black")+scale_fill_gradient2(low = "red", high = "blue")+ggtitle("theta chr 21 backward")+xlab("strand states")+ylab("cells")

grid.arrange(theta.chr1.0.plt, theta.chr1.16.plt, nrow=1, ncol=2)


## looking into counts data

plot_w_c_heatmaps <- function(cl1, cl2, chrom){

	for (counts in counts.l)
		colnames(counts) <- c('w', 'c')

	back.clust <- which(hard.cluster$unmerged.ord==cl1)
	back.clust.chr1.16 <- which(hard.cluster$unmerged.ord==cl1 & hard.cluster$pb.chr=='chr1' & hard.cluster$pb.flag==16)
	fore.clust <- which(hard.cluster$unmerged.ord==cl2)
	fore.clust.chr1.0 <- which(hard.cluster$unmerged.ord==cl2 & hard.cluster$pb.chr=='chr1' & hard.cluster$pb.flag==0)

	counts.back.clust <- lapply(counts.l, function(x) x[back.clust,])
	counts.back.clust.chr1.16 <- lapply(counts.l, function(x) x[back.clust.chr1.16,])
	counts.fore.clust <- lapply(counts.l, function(x) x[fore.clust,])
	counts.fore.clust.chr1.0 <- lapply(counts.l, function(x) x[fore.clust.chr1.0,])


	strand.states = c('ww', 'cc', 'wc')

	estim.states.back.clust <- lapply(theta.l, function(x) strand.states[which.max(x[cl1,])])
	estim.states.fore.clust <- lapply(theta.l, function(x) strand.states[which.max(x[cl2,])])


	plots = list()
	for (cell in 1:length(counts.l)){
		cell.name <- names(counts.l)[cell]
		wc.counts.back.clust <- as.matrix(table(data.table(counts.back.clust[[cell]])))
		wc.counts.back.clust[1,1] <- 0
		wc.counts.back.clust.chr1.16 <- as.matrix(table(data.table(counts.back.clust.chr1.16[[cell]])))
		wc.counts.back.clust.chr1.16[1,1] <- 0
		wc.counts.fore.clust <- as.matrix(table(data.table(counts.fore.clust[[cell]])))
		wc.counts.fore.clust[1,1] <- 0
		wc.counts.fore.clust.chr1.0 <- as.matrix(table(data.table(counts.fore.clust.chr1.0[[cell]])))
		wc.counts.fore.clust.chr1.0[1,1] <- 0

		plots[[length(plots)+1]] <- pheatmap(wc.counts.back.clust, color = colorRampPalette(c("white", "red"))(100), cluster_rows = F,cluster_cols = F, main=cell.name)[[4]]
		plots[[length(plots)+1]] <- pheatmap(wc.counts.back.clust.chr1.16, color = colorRampPalette(c("white", "red"))(100), cluster_rows = F,cluster_cols = F, main=estim.states.back.clust[[cell]])[[4]]
		plots[[length(plots)+1]] <- pheatmap(wc.counts.fore.clust, color = colorRampPalette(c("white", "red"))(100), cluster_rows = F,cluster_cols = F, main=cell.name)[[4]]
		plots[[length(plots)+1]] <- pheatmap(wc.counts.fore.clust.chr1.0, color = colorRampPalette(c("white", "red"))(100), cluster_rows = F,cluster_cols = F, main=estim.states.fore.clust[[cell]])[[4]]
	}


	pdf(paste0('counts_clust', cl1, '_clust', cl2, '.pdf'))

	for (i in 0:(ceiling(length(plots)/16)-1)){
		do.call(grid.arrange, c(plots[(i*16+1):min(length(plots),(i+1)*16)], nrow = 4, ncol = 4))
	}

	dev.off()
}

#plot_w_c_heatmaps(21, 63) # chr1
plot_w_c_heatmaps(6, 9)

# TODO: create a PDF file, put every four cells in one page... add estimated strand state to the title as well

### TODO: figure out what are the theta columns
## If third column is wc, there are a lot of cells with wc type in chr21 (more than 80%) !!!
## double check the theta estimate function
## check for the number of reads falling into repeat regions

## maybe the SS to PB mapping has a very low quality


## try to prepare some plots for showing the theta params of different clusters colored (shown) by their ground truth chrom_dirs

## If the problem arises from the cluster sizes, maybe clustering the clusters by thier theta_wc s first to the #chroms+1 works better since we double the sizes of clusters, and it makes the small clusters bigger and more visible

## also plot the theta_wc s of different clusters colored (shown) by their ground truth chroms





