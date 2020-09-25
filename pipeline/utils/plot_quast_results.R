library(data.table)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)

chr.names <- paste0('chr', c(1:22, 'X'))

quast.dir = 'aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/split/haplo_split_min_2_bubbles_k63_a3_l23/canu/quast_results'
dir.names <- list.files(quast.dir)
dir.names <- data.table(dirname=dir.names)
dir.names[, ref_haplo:=strsplit(dirname, '_')[[1]][2], by=dirname]
dir.names[, clust.forward:=strsplit(dirname, '_')[[1]][3], by=dirname]
dir.names[, clust.backward:=strsplit(dirname, '_')[[1]][4], by=dirname]
dir.names[, assembly_haplo:=strsplit(dirname, '_')[[1]][5], by=dirname]
dir.names[assembly_haplo=="h2.p", assembly_haplo:="H2"]
dir.names[assembly_haplo=="h1.p", assembly_haplo:="H1"]
dir.names[ref_haplo=="H1", ref_haplo:="H2"]
dir.names[ref_haplo=="H0", ref_haplo:="H1"]

clust.pairs <- fread('aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/clust_partners.txt')
clust.pairs[, original.chrom:=strsplit(original.chrom, '_')[[1]][1], by=original.chrom]

dir.names <- merge(dir.names, clust.pairs[, .(clust.forward, original.chrom)], by='clust.forward')

all.results <- data.table()
for (f in dir.names[, dirname]){
	results <- fread(file.path(quast.dir, f, 'report.tsv'), header=F)
	results[, dirname:=f]
	results <- results[V1 %in% c('N50', 'NA50', 'NA75', 'LA50', 'LA75')]
	results[, V2:=lapply(.SD, as.integer), .SDcols="V2"]
	results <- data.table::dcast(results, dirname ~ V1, value.var="V2")
	all.results <- rbind(all.results, results, fill=TRUE)
}

dir.names <- merge(dir.names, all.results, by='dirname', all=T)

dir.names[ref_haplo==assembly_haplo, sum_NA50:=sum(NA50), by=original.chrom]
dir.names[ref_haplo!=assembly_haplo, sum_NA50:=sum(NA50), by=original.chrom]
dir.names[, max_sum_NA50:=max(sum_NA50), by=original.chrom]

dir.names <- dir.names[sum_NA50==max_sum_NA50]
dir.names[, `:=`(sum_NA50=NULL, max_sum_NA50=NULL)]

ord.chr <- factor(dir.names$original.chrom, levels=chr.names)

na50.plt <- ggplot(dir.names, aes(x=ord.chr, y=NA50, fill=ref_haplo, width=0.5)) + geom_bar(stat="identity", color="black", position=position_dodge()) + scale_fill_hue(l=40) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
n50.plt <- ggplot(dir.names, aes(x=ord.chr, y=N50, fill=ref_haplo, width=0.5)) + geom_bar(stat="identity", color="black", position=position_dodge()) + scale_fill_hue(l=40) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
la50.plt <- ggplot(dir.names, aes(x=ord.chr, y=LA50, fill=ref_haplo, width=0.5)) + geom_bar(stat="identity", color="black", position=position_dodge()) + scale_fill_hue(l=40) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

n50 <- grid.arrange(n50.plt, na50.plt)

ggsave('assembly_results/canu/N50.pdf', grid.arrange(n50.plt, na50.plt))
ggsave('assembly_results/canu/LA50.pdf', la50.plt)


