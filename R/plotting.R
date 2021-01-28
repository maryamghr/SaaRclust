#' Plot data quality measures.
#'
#' Takes imported data table and plots some relevant data quality measures.
#'
#' @param summary.tab Imported data table.
#' @author David Porubsky
#' @export

plotQualMeasure <- function(summary.tab) {

  mapp.dist.tab <- summary.tab$mapp.stat.counts
  mapp.gaps.tab <- summary.tab$mapp.gaps.stat
  SScov.stat.m <- summary.tab$SScov.stat          
  ord <- order(match(rownames(SScov.stat.m), mapp.dist.tab$PBreadNames))
  SScov.stat.m <- SScov.stat.m[ord,]
  
  #plotting distribution of SSreads mapped to PB reads
  read.map.dist.mean <- round(mean(mapp.dist.tab$SSread.perPB))
  read.dist.plt <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSread.perPB)) + theme_bw() + geom_hline(yintercept = read.map.dist.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), read.map.dist.mean, label = read.map.dist.mean, vjust = -1), color="red") + xlab("Sorted PB reads") + ylab("# of mapped SS reads")
  read.dist.plt.log <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSread.perPB)) + theme_bw() + geom_hline(yintercept = read.map.dist.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), read.map.dist.mean, label = read.map.dist.mean, vjust = -1), color="red") + xlab("Sorted PB reads") + ylab("# of mapped SS reads (log10)") + scale_y_continuous(trans = "log10")
  
  #plotting distribution of SSlibs represented per PB read
  SSlib.perPB.mean <- round(mean(mapp.dist.tab$SSlib.perPB))
  SSlib.perPB.plt <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSlib.perPB)) + theme_bw() + geom_hline(yintercept = SSlib.perPB.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), SSlib.perPB.mean, label = SSlib.perPB.mean, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("# of SS libs per PB read")
  SSlib.perPB.plt.log <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSlib.perPB)) + theme_bw() + geom_hline(yintercept = SSlib.perPB.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), SSlib.perPB.mean, label = SSlib.perPB.mean, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("# of SS libs per PB read (log10)") + scale_y_continuous(trans = "log10")
  
  #plotting sum of gaps per PB read normalized by SSread counts mapped to a given PB read
  gaps.perPB.top1perc <- round(quantile(mapp.dist.tab$gaps.perPB.norm,prob=1-1/100)) #99th quantile
  gaps.perPB.plt <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=gaps.perPB.norm)) + theme_bw() + geom_hline(yintercept = gaps.perPB.top1perc, color="red") + geom_text(aes(nrow(mapp.dist.tab), gaps.perPB.top1perc, label = gaps.perPB.top1perc, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("Normalized sum of gaps per PB read")
  
  #plotting counts of SS reads per lib per PB read
  mapp.SScov.stat.mean <- mean(SScov.stat.m[SScov.stat.m>0])
  SScov.stat.m[SScov.stat.m==0] <- NA
  SScov.stat.df <- data.frame(value=rowMeans(SScov.stat.m, na.rm = T))
  SScov.stat.plt <- ggplot(SScov.stat.df) + geom_line(aes(x=c(1:nrow(SScov.stat.df)),y=value)) + theme_bw() + geom_hline(yintercept = mapp.SScov.stat.mean, color="red") + geom_text(aes(nrow(SScov.stat.df), mapp.SScov.stat.mean, label = mapp.SScov.stat.mean, vjust = -1), color="red") + xlab("PB reads sorted by the mean number\nof SSreads per SSlib per PBread") + ylab("Mean counts of SSreads per SSlib per PB read")
  
  #plotting match with gaps per SS read flagged TRUE or FALSE based on mapping agreement
  gaps.perSS.mean <- round(mean(mapp.gaps.tab$matchWithgaps))
  accur.counts <- table(mapp.gaps.tab$mapp.accur)
  falses <- round((accur.counts[1]/sum(accur.counts))*100, digits = 2)
  trues <- round((accur.counts[2]/sum(accur.counts))*100, digits = 2)
  false.lab <- paste('MISS ', falses, '%', sep = "")
  true.lab <- paste('MATCH ', trues, '%', sep = "")
  gaps.perSS.plt <- ggplot(mapp.gaps.tab) + geom_linerange(aes(x=c(1:nrow(mapp.gaps.tab)),ymin=0, ymax=matchWithgaps, color=mapp.accur)) + theme_bw() + geom_hline(yintercept = gaps.perSS.mean, color="red") + geom_text(aes(nrow(mapp.gaps.tab), gaps.perSS.mean, label = gaps.perSS.mean, vjust = -1), color="red") + xlab("Sorted SS reads by the length\nof SS read alignment (bp)") + ylab("(bp) SS read alignment with gaps (log10)") + scale_y_continuous(trans = "log10") + scale_color_manual(labels = c(false.lab, true.lab), values = c("red", "green"))
 
  #plot PB read length distribution
  PBreadLenDist <- ggplot(summary.tab$PBreadLenDist) + geom_linerange(aes(x=midpoints, ymin=0, ymax=freq), size=3)
  
  suppressWarnings( plt <- plot_grid(read.dist.plt, read.dist.plt.log, SSlib.perPB.plt, SSlib.perPB.plt.log, SScov.stat.plt, gaps.perPB.plt, PBreadLenDist, ncol = 2) )
  return(plt)
}


#' Plot heatmap of responsibilities of each PB reads for each cluster as a probability value.
#'
#' @param pVal.df A \code{data.frame} of probability values.
#' @param colOrder A \code{vector} of indices representing of cluster order on heatmap.
#' @param num.clusters Number of cluster present in probability table.
#' @author David Porubsky
#' @export

plotHeatmap <- function(soft.clust=NULL, colOrder=NULL, num.clusters=NULL) {

  #order clusters based on most likely chromosome partners                               
  if (!is.null(colOrder)) {                               
    soft.clust$soft.pVal[,1:length(colOrder)] <- soft.clust$soft.pVal[,colOrder]
    colnames(soft.clust)[1:num.clusters] <- colOrder
  }
  
  if (!is.null(num.clusters)) {

    chr.ids <- names(sort(table(soft.clust$chrom), decreasing = T))
    chr.ids <- gsub('^chr', '', chr.ids)
    chr.ids <- sort(as.numeric(chr.ids))
    chr.colors <- rep(c("gray48","gray72"), ceiling(length(chr.ids)/2))
    chr.colors <- chr.colors[1:length(chr.ids)]
    names(chr.colors) <- chr.ids
    
    soft.clust$chrom <- gsub('^chr', '', soft.clust$chrom)
    soft.clust$chrom <- factor(soft.clust$chrom, levels=chr.ids)
    soft.clust$soft.pVal <- soft.clust$soft.pVal[order(soft.clust$chrom),]
    # soft.clust$chrom <- soft.clust$chrom[order(soft.clust$chrom)]
    
    #set unexpected directionality flags to 1
    soft.clust$flag[soft.clust$flag != 0 & soft.clust$flag != 16] <- 1
    
    ha1 = rowAnnotation(df = data.frame(chr = soft.clust$chrom), col = list(chr=chr.colors))
    ha2 = rowAnnotation(df = data.frame(dir = soft.clust$flag), col = list(dir = c('0'="chocolate1", '16'="chartreuse3", '1'="white")))
    hm <- Heatmap(soft.clust$soft.pVal[,c(1:num.clusters)], name = "Probs", cluster_columns = F, cluster_rows = F, show_row_names = FALSE)
    #hm <- Heatmap(soft.clust$soft.pVal[,c(1:num.clusters)], name = "Probs", cluster_columns = F, cluster_rows = F, show_row_names = FALSE)
    hm + ha1 + ha2
    return(hm)
  } else {
    message("num.clusters not specified!!!\n")
  }  
}  


#' Plot accuracy of clustring in respect to known location of each PacBio read.
#'
#' @inheritParams plotHeatmap
#' @param thresh A \code{vector} of probability thresholds to obtain accuracy of clustering. 
#' @author David Porubsky
#' @export

plotClustAccuracy <- function(pVal.df=NULL, chrom=NULL, flag=NULL, num.clusters=NULL, thresh=c(0.5,0.6,0.7,0.8,0.9,0.99)) {
  pVals <- pVal.df[,c(1:num.clusters)]
  chr.rows <- chrom
  
  acc.l <- list()
  for (th in thresh) {
    max.pVal <- apply(pVals, 1, max)
    mask <- max.pVal >= th
    ord <- apply(pVals[mask,], 1, which.max)
    chr.clusts <- split(chr.rows[mask], ord)
    chr.clusts <- lapply(chr.clusts, unlist) #Check why is this problem???
    clust.acc <- getClusterAcc(chr.clusts)
    frac.corr <- sum(clust.acc$stat$trues)/(sum(clust.acc$stat$trues) + sum(clust.acc$stat$miss))
    eval.reads <- table(mask)
    eval.reads <- eval.reads[2]/sum(eval.reads)
    acc.l[[as.character(th)]] <- c(frac.corr, eval.reads)
  } 
  
  #get hard clust accuracy
  chr.clusts <- split(chr.rows, pVal.df$hardClust)
  chr.clusts <- lapply(chr.clusts, unlist) #Check why is this problem???
  clust.acc <- getClusterAcc(chr.clusts)
  frac.corr <- sum(clust.acc$stat$trues)/(sum(clust.acc$stat$trues) + sum(clust.acc$stat$miss))
  
  m <- do.call(rbind, acc.l)
  df <- data.frame(values=m[,1], eval=m[,2], thresh = factor(thresh, levels=rev(thresh)))
  HC <- data.frame(values=frac.corr, eval=1, thresh='HardClust')
  
  plt <- ggplot(df) + geom_point(aes(x=values, y=eval, color=thresh), size=5) + geom_linerange(aes(ymin=-Inf, x=values, ymax=eval, color=thresh)) + scale_y_continuous(limits = c(0,1)) + scale_color_manual(values = brewer.pal(n=7, name="Set1"), name="Prob threshold") + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads")
  plt <- plt + geom_point(data=HC, aes(x=values, y=eval, color=thresh), size=5, inherit.aes = F) + geom_linerange(data=HC, aes(ymin=-Inf, x=values, ymax=eval, color=thresh), inherit.aes = F)
  return(plt)
}


#' Plot theta estimates resulting from EM algorithm.
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @param title A \code{character} to use as a title of the plot.
#' @importFrom reshape2 melt
#' @author David Porubsky
#' @export

plotThetaEstimates <- function(theta.param=NULL, title=NULL) {
  
  plt.data <- list()
  for (j in 1:length(theta.param)) {
    df <- as.data.frame(theta.param[[j]])
    df$clustID <- rownames(df)
    df.plt <- suppressMessages( reshape2::melt(df) )
    df.plt$cell <- j
    plt.data[[j]] <- df.plt
  }
  plt.data.df <- do.call(rbind, plt.data)

  my_theme <-  theme(panel.spacing = unit(0, "lines"), 
                   strip.text.y = element_text(angle = 0),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())
  if (is.null(title)) {
    plt <- ggplot(plt.data.df , aes(x=clustID, y=value, fill=variable)) + geom_bar(stat='identity', width=1) + facet_grid(cell ~ .) + scale_fill_manual(values = c('prob.cc'="paleturquoise4", 'prob.mix'="olivedrab",'prob.ww'="sandybrown")) + my_theme
  } else {
    plt <- ggplot(plt.data.df , aes(x=clustID, y=value, fill=variable)) + geom_bar(stat='identity', width=1) + facet_grid(cell ~ .) + scale_fill_manual(values = c('prob.cc'="paleturquoise4", 'prob.mix'="olivedrab",'prob.ww'="sandybrown")) + ggtitle(title) + my_theme
  }
  return(plt)
}


#' Plot distribution of short reads mapped on top of PB reads
#'
#' @param count.list A \code{list} of short read mappings per library.
#' @author David Porubsky
#' @export

plotReadMappingDist <- function(count.list=NULL) {
  
  SSperPB <- list()
  for (j in 1:length(count.list)) {

      lib.aligns <- count.list[[j]]
      counts <- table(lib.aligns$PBreadNames)
      SSperPB[[j]] <- counts
  }
  all.counts <- do.call(rbind, SSperPB)
  plt.df1 <- as.data.frame(table(all.counts))
  plt1 <- ggplot(plt.df1) + geom_bar(aes(x=as.numeric(all.counts), y=Freq), stat='identity', fill="red") + xlab("# of ShortReads per PBread per Library") + ylab("Frequency") + scale_x_continuous(breaks = as.numeric(plt.df1$all.counts), labels = plt.df1$all.counts)
  
  count.list.collapsed <- do.call(rbind, count.list)
  counts <- table(count.list.collapsed$PBreadNames)
  plt.df2 <- as.data.frame(table(counts))
  
  is.odd <- function(x) x %% 2 != 0
  breaks <- as.numeric(plt.df2$counts)[ is.odd(as.numeric(plt.df2$counts)) ]
  plt2 <- ggplot(plt.df2) + geom_bar(aes(x=as.numeric(counts), y=Freq), stat='identity', fill="red") + xlab("# of ShortReads per PBread") + ylab("Frequency") + scale_x_continuous(breaks = breaks, labels = breaks)
  
  plt <- plot_grid(plt1, plt2, nrow = 1, rel_widths = c(1,2))
  return(plt)
}


#' Plot coverage of short reads mapped on top of PB reads
#'
#' @param minimap.tab A \code{data.frame} of short read mappings per PacBio read in maf.
#' @author David Porubsky
#' @export

plotReadAlignments <- function(minimap.tab=NULL) {
  #Convert table of alignments into GRanges object and then split into GRangesList by StrandS library ID
  minimap.tab.gr <- GenomicRanges::GRanges(seqnames=minimap.tab$PBchrom, strand=minimap.tab$strand, ranges=IRanges(start=minimap.tab$TargetCoordStart, end=minimap.tab$TargetCoordend), PBreadLen=minimap.tab$PBreadLen, SSlibNames=minimap.tab$SSlibNames)
  minimap.tab.grl <- GenomicRanges::split(minimap.tab.gr, minimap.tab.gr$SSlibNames)
  
  #get the name of PB read
  readID <- as.character(unique(minimap.tab$PBreadNames))
  
  all.libs <- list()
  #probs.l <- list()
  for (i in 1:length(minimap.tab.grl)) {
    gr <- minimap.tab.grl[[i]]
    gr$level <- GenomicRanges::disjointBins(gr)
    gr$level[which(GenomicRanges::strand(gr) == '-')] <- gr$level[which(GenomicRanges::strand(gr) == '-')] * -1
    
    #Get probabilities for StrandS read distribution
    dirRead.counts <- table(GenomicRanges::strand(gr))
    probs <- countProb(minusCounts = dirRead.counts['-'], plusCounts = dirRead.counts["+"], alpha = 0.1)
    probs.norm <- probs/sum(probs) #normalize prob values to 1
    probs.string <- paste(probs.norm, collapse = ", ")
    gr$probs <- probs.string
    
    #probs.df <- data.frame(minus=dirRead.counts['-'], plus=dirRead.counts["+"], ww=probs[,1], cc=probs[,2] ,wc=probs[,3], max=which.max(probs))
    #probs.l[[i]] <- probs.df
    
    plt.df <- as.data.frame(gr)
    all.libs[[i]] <- plt.df
  }
  all.libs.df <- do.call(rbind, all.libs)
  #all.probs.df <- do.call(rbind, probs.l)
  
  readLen <- data.frame(start=0, end=unique(all.libs.df$PBreadLen))
  plt <- ggplot(all.libs.df) + geom_linerange(data=readLen, aes(x=0, ymin=start, ymax=end), color="black") + geom_linerange(aes(x=level, ymin=start, ymax=end, color=strand)) + coord_flip() + scale_color_manual(values = c("paleturquoise4","sandybrown")) + xlab("") + facet_grid(SSlibNames ~ ., scales = 'free') + geom_text(aes(x=Inf,y=0, vjust=1, hjust=0), label=all.libs.df$probs) + ggtitle(readID) + theme(strip.text.y = element_text(angle = 360))
  return(plt)
}

#' Plots heatmap of w reads fractions in chrom/libs
#'
#' @param counts.l A list of W/C read count data tables (rows=long reads, cols=w,c) per single-cell
#' @author Maryam Ghareghani
#' @export


plot_heatmap_wfrac <- function(counts.l, compute_Wfrac=F, by_chrom=F)
{
  # avg.unitigs.w.frac is from count.selected obtained from get_representative_alignments in import.R
  
  w.frac <- counts.l
  
  if (compute_Wfrac)
  {
    w.frac[, W.frac:=(w-c)/(w+c)]
    #w.frac[, `:=`(w=NULL, c=NULL)]
  }
  
  mat <- dcast(w.frac, rname+chrom+flag~lib, value.var='W.frac')
  # test
  #mat <- mat[1:100,1:10]
  
  chroms <- paste0('chr',c(1:22, 'X'))
  if (by_chrom) {
    expected.clusters <- chroms
  } else {
    expected.clusters <- paste0(rep(chroms, each=2), rep(c('_0','_16'),23))
    mat[, chrom:=paste0(chrom,'_',flag)]
    
  }
  row.order <- factor(expected.clusters, levels = expected.clusters)

  mat <- mat[order(factor(chrom, levels=row.order))]
  
  annot <- HeatmapAnnotation(data.frame(chromosome=mat[, chrom]), which = 'row')
  if (!by_chrom){
    annot <- annot + HeatmapAnnotation(data.frame(flag=as.character(mat[, flag])), which = 'row')
  }
  
  row.names <- mat[, rname]
  mat[, `:=`(rname=NULL,chrom=NULL,flag=NULL)]
  # only cluster by chrom:
  #  mat <- abs(mat)
  mat <- as.matrix(mat)
  rownames(mat) <- row.names
  
  if (by_chrom) {mat <- abs(mat)}
  
  d= as.matrix(dist(mat))
  
  plt.clust <- Heatmap(mat, name='w-c fraction', show_row_names = F, show_column_names = F) + annot #+ chrom.annot + flag.annot
  plt.clust.rows <- Heatmap(mat, name='w-c fraction', cluster_columns = F, show_row_names = F, show_column_names = F) + annot #+ chrom.annot + flag.annot
  plt <- Heatmap(mat, name='w-c fraction', cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F) + annot #+ chrom.annot + flag.annot
  plt.dist.clust <- Heatmap(d, name='dist of w-c fraction vectors', show_row_names = F, show_column_names = F) + annot #+ chrom.annot + flag.annot
  plt.dist <- Heatmap(d, name='dist of w-c fraction vectors', cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F) + annot
  
  list(plt.clust=plt.clust, plt.clust.rows=plt.clust.rows, plt=plt, plt.dist.clust=plt.dist.clust, plt.dist=plt.dist, features.matrix=mat)
}

# hard clust heatmap

clust.chrom.count <- data.table(chrom=hard.clust$chrom, clust=hard.clust$ord)[!is.na(chrom)]
clust.chrom.count[, count:=.N, by=.(chrom, clust)]
clust.chrom.count <- clust.chrom.count[, head(.SD, 1), by=.(chrom, clust)]
mat <- dcast(clust.chrom.count, chrom~clust, value.var = 'count', fill=0)
mat <- mat[order(factor(chrom, levels=paste0('chr', c(1:22, 'X'))))]
row.names <- mat[, chrom]
mat[, chrom:=NULL]
mat <- as.matrix(mat)
rownames(mat) <- row.names
heatmap(mat, scale = 'column', Rowv = NA, Colv = NA, col= colorRampPalette(brewer.pal(8, "Blues"))(25))

# input data insights- bar plot of true clusters in all and selected unitigs
ggplot(chrom.flag)+ geom_bar(aes())

# complex Heatmap:
#Heatmap(mat)
################### SS Wfrac in long reads/unitigs
counts <- data.table()
for (i in 1:length(counts.l)){
  counts.dt <- ratios.l[[i]] # counts.l[[i]]
  lib.name=names(counts.l)[i]
  rnames=rownames(counts.l[[i]]) #(counts.dt)
  
  counts.dt[, `:=`(rname=rnames, lib.name=lib.name)]
  
  counts <- rbind(counts, counts.dt)
}

colnames(counts)[1] <- 'W.frac' #counts[, W.frac:=(w-c)/(w+c)]
counts.selected <- counts

rnames <- counts[, unique(rname)]
lib.names <- counts[, unique(lib.name)]
expand.counts <- data.table(expand.grid(rnames, lib.names))
colnames(expand.counts) <- c('rname', 'lib.name')
counts <- merge(counts, expand.counts, by=c('rname', 'lib.name'), all=T)
counts[is.na(counts)] <- 0

counts <- merge(counts, chrom.flag, by='rname')

# subsetting the highest coverage unitigs
counts[, ss.cov:=sum(w+c), by=rname]
setkey(counts, ss.cov)

counts.selected <- counts[, head(.SD, 1), by=rname]
counts.selected <- counts.selected[ss.cov >= min.ss.cov]

selected.rname <- counts.selected[, rname]
if (nrow(counts.selected) > num.alignments){
  selected.rname <- tail(selected.rname, num.alignments)
}
counts.selected <- counts[rname %in% selected.rname]
counts.selected[, avg.w.frac:=mean(W.frac[W.frac!=0]), by=.(lib.name,chrom,flag)]

avg.unitigs.w.frac <- counts.selected[, head(.SD,1), by=.(lib.name,chrom,flag)][, .(lib=lib.name, chrom, flag, W.frac=avg.w.frac)]
avg.unitigs.w.frac[, chrom:=paste0(chrom,'_', flag, '_unitig'), by=.(chrom, flag)]
avg.unitigs.w.frac[, flag:=NULL]

w.frac = rbind()