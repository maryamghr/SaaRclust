setwd("/home/maryam/research/haploclust/haploclust/pipeline/")
.libPaths(c(.libPaths(), 'utils/R-packages/'))
library(SaaRclust)
library(data.table)
library(Rsamtools)
library(doParallel)
library(matrixStats)
library(ComplexHeatmap)
source('../R/timedMessage.R')
source('../R/import.R')
source('../R/hardClust.R')
source('../R/SaaRclust.R')
source('../R/wrapper.R')
source('../R/EMclust.R')
source('../R/calcProbs.R')
source('../R/helperFuctions.R')
source('utils/evaluate_theta_param.R')
# contibait: http://bioconductor.org/packages/release/bioc/manuals/contiBAIT/man/contiBAIT.pdf
inputfolder <- '../../HG00733/bwa_ss_unitigs'
outputfolder <- '../../HG00733/SaaRclust'
input_type <- 'bam'
num.clusters <- 80
EM.iter <- 2
alpha <- 0.01
minLib <- 35
upperQ <- 1
logL.th <- 1
HC.only = FALSE
store.counts = FALSE
store.bestAlign = TRUE
numAlignments = 30000
log.scale=TRUE
numCPU=32
fileID='HG00733'
args <- c(inputfolder, outputfolder, input_type, num.clusters, alpha, numAlignments, log.scale, EM.iter, minLib, upperQ, logL.th, fileID, numCPU)
cellNum <- NULL
theta.constrain=FALSE

outputfilename=paste0("hardClusteringResults_",num.clusters,"clusters.RData")

clust <- runSaaRclust(inputfolder = args[1], outputfolder = args[2], input_type=args[3], num.clusters=as.numeric(args[4]),
                           EM.iter=args[8], alpha=as.numeric(args[5]), minLib=as.numeric(args[9]), upperQ=as.numeric(args[10]),
                           logL.th=args[11], theta.constrain=FALSE, store.counts=FALSE, store.bestAlign=TRUE, numAlignments=as.numeric(args[6]),
                           HC.only=FALSE, log.scale=args[7], outputfilename=outputfilename, fileID=args[12], numCPU=as.numeric(args[13]))

hard.clust <- clust[[1]]
soft.clust <- clust[[2]]

utg.bam.file <- '../../HG00733/hifiasm/ref_aln/asm.r_utg.bam'

hard.clust <- addChromFlag(hard.clust, utg.bam.file, 'hard')
soft.clust <- addChromFlag(soft.clust, utg.bam.file, 'soft')

# numfoundclusters for hard.clust...
# check random shuffling
# check hard clustering for all reads

 
# counting w/c in bam files
# inputfolder <- '../../HG00733/bwa_ss_unitigs/'
bam.files <- list.files(path=inputfolder, pattern='.mdup.bam$', full.names=TRUE)
counts.l <- count.wc.bam(bam.files, numCPU = 42)
counts.l <- get_representative_counts(counts.l, numAlignments)

### Perform k-means hard clustering method ###
set.seed(1000) #in order to reproduce hard clustering results
hardClust.ord <- hardClust(counts.l, num.clusters=num.clusters, num.iter = 10)
names(hardClust.ord) <- rownames(counts.l[[1]])

# looking into hard clustering features
num.cells <- ncol(ratios.m)
colnames(ratios.m) <- paste0('cell', 1:num.cells)
ratios.m[, rname:=rownames(counts.l[[1]])]
ratios.m <- merge(ratios.m, chrom.flag, by='rname')
ratios.m[, chrom_flag:=paste0(chrom, '_',  flag), by=.(chrom, flag)]

theta.estim <- estimateTheta(counts.l, hard.clust=hardClust.ord, alpha=alpha)
hardClust.ord.merged <- mergeClusters(hard.clust=hardClust.ord, theta.l=theta.estim, k=47)
#Re-estimate theta parameter after cluster merging
theta.estim <- estimateTheta(counts.l, hard.clust=hardClust.ord.merged, alpha=alpha)

#Initialize theta parameter
theta.param <- theta.estim
#Estimate pi parameter based on # of PB reads in each cluster
readsPerCluts <- table(hardClust.ord.merged)
pi.param <- readsPerCluts/sum(readsPerCluts)

#save hard clustering results into a file
hard.clust <- list(ord=hardClust.ord.merged, unmerged.ord=hardClust.ord, theta.param=theta.param)#, pi.param=pi.param, pb.chr = chr.rows, pb.flag = pb.flag)
destination <- file.path(Clusters.store, outputfilename)

################## eval:
hard.clust <- get(load('../../HG00733/SaaRclust/Clusters/hardclusters.RData'))
soft.clust <- get(load('../../HG00733/SaaRclust/Clusters/soft_clusters.RData'))
bam.file <- '../../HG00733/hifiasm/ref_aln/asm.r_utg.haplotagged.bam'

#chrom.flag <- fread('../../HG00733/SaaRclust/chrom_flags.data')
counts.l <- get(load('../../HG00733/SaaRclust/read_counts.RData'))
chrom.flag <- getChromFlag(bam.file)
hard.clust <- addChromFlag(clust.obj=hard.clust, bam.file=bam.file, clust.type = 'hard')
soft.clust <- addChromFlag(soft.clust, bam.file=bam.file, clust.type = 'soft')
valid.chroms <- paste0('chr', c(1:22, "X"))
expected.clusters <- paste0(rep(valid.chroms, each=2), rep(c('_0','_16'),23))


numFoundClusters(hard.clust$ord, hard.clust$chr, hard.clust$flag, expected.clusters) 
# hard clust
bam.files <- list.files(path=inputfolder, pattern='.mdup.bam$', full.names=TRUE)
#counts.l <- count.wc.bam(bam.files, numCPU = 42)
counts.representative.1000 <- get_representative_counts(counts.l=counts.l, num.alignments = 1000, min.ss.cov=200, plot.hist=TRUE, chrom.flag=chrom.flag)
counts.representative <- get_representative_counts(counts.l=counts.l, num.alignments = 10000, min.ss.cov=200, plot.hist=TRUE, chrom.flag=chrom.flag)
counts.selected.l <- counts.representative[[2]]
counts.selected <- counts.representative[[3]]

#purifying counts.selected
# counts.selected[, non.zero.cov:=0][w+c>0, non.zero.cov:=1]
# counts.selected[, non.zero.cov:=sum(non.zero.cov), by=rname]
# counts.selected[, `:=`(n.ww=0, n.cc=0)]
# counts.selected[w/(w+c) < 0.2, n.cc:=1]
# counts.selected[w/(w+c) > 0.8, n.ww:=1]
# counts.selected[, n.ww.cc:=sum(n.ww+n.cc), by=rname]
# filt <- counts.selected[, head(.SD), by=rname]
# filt <- merge(filt, chrom.flag, by='rname', all.x=T)
# filt[, norm.n.ww.cc:=n.ww.cc/non.zero.cov]
# 
# ggplot(filt)+geom_histogram(aes(norm.n.ww.cc, fill=is.na(chrom)), binwidth=0.01)
# ggplot(filt[norm.n.ww.cc>0.25 & norm.n.ww.cc<0.75])+geom_bar(aes(x=paste0(chrom, '_', flag)))
# 
# filt.rname <- filt[norm.n.ww.cc>0.25 & norm.n.ww.cc<0.75, rname]
# 
# counts.selected <- counts.selected[rname %in% filt.rname, .(rname,w,c,lib.name, ss.cov)]


### hclust
wfrac.matrix <- get.wfrac.matrix(counts.dt=counts.selected)
hc <- hard.hclust(wfrac.matrix=wfrac.matrix, num.clusters=60)

## plot heatmap
counts.selected <- merge(counts.selected, chrom.flag, by='rname', all.x=T)
plt.l <- plot_heatmap_wfrac(input.dt=counts.selected, compute_Wfrac=F)
plt.chrom.l <- plot_heatmap_wfrac(input.dt=counts.selected, by_chrom =F)

# num.clusters=47: missing clusters: chr14_0, chr22_0, clust accuracy: 0.9745495
# num.clusters=50: missing clusters: chr22_0, clust accuracy: 0.9877252
# num.clusters=60,80 on 10000 reads yields very good results: 99.8% acc with finding all chrom_dir clsuters, time=18s

hc.chrom <- hard.hclust(wfrac.matrix=wfrac.matrix, num.clusters=30, by_chrom = T) # with 30 clusters on 10000 reads finds all chromosomes, acc=99.6%
clust <- data.table(rname=names(hc), clust=hc)
clust <- merge(clust, chrom.flag, by='rname', all.x=T)
clust.by.chrom <- data.table(rname=names(hc.chrom), clust=hc.chrom)
clust.by.chrom <- merge(clust.by.chrom, chrom.flag, by='rname', all.x=T)

clust.to.chrom <- numFoundClusters(ord=clust$clust, chr=clust$chrom, flag=clust$flag, expected.clusters)
numFoundChromosomes(ord=clust.by.chrom$clust, clust$chrom, valid.chroms)

clust[, num.clust:=.N, by=clust]
# computing cluster counts
clust.num <- clust[, .(rname,clust,num.clust)]
clust.num <- merge(clust.num, clust.to.chrom, by='clust')
ggplot(clust.num) + geom_bar(aes(x=clust, fill=is.na(chrom)))
clust.num <- clust.num[, head(.SD,1), by=clust]
setkey(clust.num, chrom, flag)
clust.num

ord <- hc

# filtering out extremly small clusters
min.clust.num <- nrow(clust)*0.002
clust <- clust[num.clust > min.clust.num]

#theta.estim <- estimateTheta(counts.selected.l, hard.clust=ord, alpha=alpha)

##### assign theta based on median wfrac values
theta <- merge(counts.selected[, .(rname,lib, W.frac)], clust[, .(rname, clust)], by='rname')
theta <- merge(theta, clust.to.chrom, by='clust', all.x=T)
theta[!is.nan(W.frac), median.w.frac:=median(W.frac), by=.(clust, lib)]
theta[, state:=0]
theta[median.w.frac <= -0.5, state:=-1, by=.(clust,lib)]
theta[median.w.frac > -0.5 & median.w.frac < 0.5, state:=0, by=.(clust,lib)]
theta[median.w.frac >= 0.5, state:=1, by=.(clust,lib)]


#theta[, `:=`(theta.ww=0, theta.cc=0, theta.wc=0)]
#theta[median.w.frac < -0.5, `:=`(theta.ww=0.05, theta.cc=0.95, theta.wc=0.05), by=.(clust,lib)]
#theta[median.w.frac > -0.5 & median.w.frac < 0.5, `:=`(theta.ww=0.05, theta.cc=0.05, theta.wc=0.95), by=.(clust,lib)]
#theta[median.w.frac > 0.5, `:=`(theta.ww=0.95, theta.cc=0.05, theta.wc=0.05), by=.(clust,lib)]

theta <- theta[, head(.SD,1), by=.(clust, lib)]
setkey(theta, clust, lib)

theta[, wc.frac:=length(which(state==0))/.N, by=clust]

wc.frac <- theta[, head(.SD, 1), by=clust]
setkey(wc.frac, chrom, flag)
wc.frac

theta.matrix <- dcast(theta, clust~lib, value.var = "state")
rownames(theta.matrix) <- theta.matrix[, clust]
theta.matrix[, clust:=NULL]
Heatmap(theta.matrix, show_row_names = T, show_column_names = F)

theta.l <- split(theta[, .(lib, state)], by='lib')
lapply(theta.l, function(dt) dt[, lib:=NULL])

# better than the probabilistic approach--- test it
# TODO: try removong the clusters with high variance

theta.estim <- split(theta[, .(lib, theta.ww, theta.cc, theta.wc)], by='lib') %>% invisible()
lapply(theta.estim, function(dt) dt[, lib:=NULL])
##

set.seed(1000)
num.clusters <- 60
km = hardClust(counts.selected.l, num.clusters= num.clusters) #, iter.max = 100)
km.chrom = hardClust(counts.selected.l, num.clusters= 40, by_chrom=T)
clust <- data.table(rname=names(km), clust=km)
clust <- merge(clust, chrom.flag, by='rname', all.x=T)
clust.by.chrom <- data.table(rname=names(km.chrom), clust=km.chrom)
clust.by.chrom <- merge(clust.by.chrom, chrom.flag, by='rname', all.x=T)
# num.clusters=50: 7 missing clusters: chr10_0, chr13_0, chr15_16, chr16_0, chr18_0, chr21_0, chr21_16, clust accuracy: 0.9027575
# num.clusters=60: missing clusters: chr21_16, chr22_16, clust accuracy: 0.9834553
# num.clusters=80 on 10000 reads yields very good results: 99.5% acc with finding all chrom_dir clsuters, time=18s

clust.to.chrom <- numFoundClusters(ord=clust$clust, chr=clust$chrom, flag=clust$flag, expected.clusters)
numFoundChromosomes(ord=clust.by.chrom$clust, clust$chrom, valid.chroms)

theta.estim <- estimateTheta(counts.selected.l, hard.clust=ord, alpha=alpha) # it has warnings

# test clust=1 (chr4_0)
test <- counts.selected[rname %in% clust[ord==1, rname]]


# theta.wc <- Reduce(c, lapply(theta.estim, function(x) x[,3]))
# table(theta.wc)
# test estimateTheta...
# problem: theta_wc is often the most likely state
# Looking at single-cell counts for cluster 1 (chr17_16) to see strand states of single-cells...
# split(counts.selected.l[[j]], hard.clust)[[1]] for j=1:115 ....
# The data (strand states) does not look clean although the cluster is rather clean containing mainly the true reads!
# double check the unitigs and ss reads mapping orientations (in IGV)
lib1.clust3 <- counts.selected.l[["HG00733.P00IL002"]][which(ord==3)] # only chr6_0
lib1.clust3[, rname:=rownames(counts.selected.l[["HG00733.P00IL002"]])[which(ord==3)]]
chrom.flag[rname %in% lib1.clust3[, rname]] # shows a mixture of chromosomes!!!
#lib1.clust3[, rname] != names(ord)[which(ord==3)]!!!

hardClust.ord.merged <- mergeClusters(hard.clust=ord, theta.l=theta.estim, k=47)
clust.to.chrom <- numFoundClusters(ord=hardClust.ord.merged, chr=clust$chrom, flag=clust$flag, expected.clusters)

theta.estim <- estimateTheta(counts.selected.l, clust.to.chrom, hard.clust=hardClust.ord.merged, alpha=alpha)
compare_theta_with_wfrac(theta=theta.estim, clust.to.chrom, countfiles.dir)
### Note: almost all thetas are wc... I guess kmeans does not work well with sparse data
# cluster only by chrom
abs.ratios <- abs(ratios.m)
abs.ratios[abs.ratios==0]=0.5
chrom.km = hardClust(counts.selected.l, num.clusters= 30)
clust <- data.table(rname=names(chrom.km), ord=chrom.km)
clust <- merge(clust, chrom.flag, by='rname', all.x=T)
#clust <- clust[chrom %in% valid.chroms]
ord = chrom.km

numFoundChromosomes(ord, chr=clust$chrom, expected.clusters = valid.chroms)
# cluster by chrom


chr.hists <- list()

for (data.file in list(hard.clust, soft.clust)){
  if ("ord" %in% names(data.file)) {utg.names <- names(data.file$ord)} else{utg.names <- rownames(data.file$soft.pVal)}
  gr.truth <- data.table(name=utg.names, chrom=data.file$chrom, flag=data.file$flag)
  gr.truth <- gr.truth[chrom %in% valid.chroms]
  gr.truth[, chrom_flag:=paste0(chrom, '_', flag), by=.(chrom, flag)]
  # order
  gr.truth[, chr.ord:=gsub('chr','',chrom), by=chrom]
  gr.truth[chrom=='chrX', chr.ord:="23"]
  gr.truth[, chr.ord:=as.numeric(chr.ord)]
  setkey(gr.truth, chr.ord)
  #######
  chr.hists[[1+length(chr.hists)]] <- ggplot(gr.truth)+geom_bar(aes(x=chrom_flag))
}

chr.flag.ord <- paste0(rep(valid.chroms, each=2), rep(c('_0','_16'),23))

## contibait
library(contiBAIT)
strand.freq <- StrandFreqMatrix(counts=as.matrix(ratios.m))
strand.states <- preprocessStrandTable(exampleStrandFreq, lowQualThreshold=0.8)
cl <- clusterContigs(strand.states$strandMatrix)
data("exampleDividedChr")
barplotLinkageGroupCalls(cl, exampleDividedChr)

# todo: PCA dim reduction:
#https://www.youtube.com/watch?v=K5HcVuLra38
pca.ratios <- prcomp(ratios.m)
summary(pca.ratios)
pca <- as.data.table(pca.ratios$x)
pca[,rname:=rownames(counts.l[[1]])]
pca = merge(pca, chrom.flag, by='rname')
pca[,chrom_flag:=paste0(chrom, '_', flag)]
ggplot(pca, aes(x=PC1, y=PC2, color=paste(chrom, flag)))+geom_point()


# plotting binned w/c read counts for large unitigs
unitigs.len <- ...
# importing bam as Granges
counts.utg000445l <- importBams(bamfolder=inputfolder, chromosomes='utg000445l', bin.length=10000)

