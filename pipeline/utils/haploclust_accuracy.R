library(data.table)
library(dplyr)
library(seqinr)
library(Rsamtools)
library(ggplot2)
library(ggvenn)
library(scales)

#setwd('/home/maryam/research/haploclust/haploclust/pipeline/')
bam.file <- "../../HG00733/hifiasm/ref_aln/asm.r_utg.haplotagged.bam"
predicted.haplo.files <- list.files('../../HG00733/phased_unitigs/', full.names = T)
predicted.chrom.file <- '../../HG00733/SaaRclust/Clusters/MLclust.data'
valid.chroms <- paste0('chr', c(1:22,'X'))

unitigs.bam <- scanBam(BamFile(bam.file), param=ScanBamParam(what=c("rname", "qname"), tag="HP", flag=scanBamFlag(isSupplementaryAlignment=FALSE)))
original.haplotypes <- data.table(ref_name=as.character(unitigs.bam[[1]]$rname),
                                  unitig_name=as.character(unitigs.bam[[1]]$qname),
                                  unitig_gr_haplo=as.character(unlist(unitigs.bam[[1]]$tag)))
original.haplotypes[!is.na(unitig_gr_haplo), unitig_gr_haplo:=paste0('H',unitig_gr_haplo)]
predicted.haplotypes <- lapply(predicted.haplo.files, fread)
predicted.haplotypes <- Reduce(rbind, predicted.haplotypes)
names(predicted.haplotypes)[1] <- 'unitig_name'
predicted.chroms <- fread(predicted.chrom.file)
predicted.chroms[, chrom_clust:=paste0('chr', chrom_clust)]
predicted.chroms[chrom_clust=='chr23', chrom_clust:='chrX']
names(predicted.chroms)[1] <- 'unitig_name'

compare.haplotypes <- merge(original.haplotypes, predicted.chroms, by='unitig_name', all.x=T)
compare.haplotypes <- compare.haplotypes[ref_name %in% valid.chroms]
# SaaRclust accuracy
message('SaaRclust Accuracy =', compare.haplotypes[ref_name==chrom_clust, .N]/compare.haplotypes[, .N])

compare.haplotypes <- merge(compare.haplotypes, predicted.haplotypes, by='unitig_name', all.x=T)

# Venn diagram for the set of haplotagged and haplo-clustered unitigs
unitigs.set <- original.haplotypes[, unitig_name]
haplotagged.unitigs <- unitigs.set %in% compare.haplotypes[!is.na(unitig_gr_haplo), unitig_name]
haploclust.unitigs <- unitigs.set %in% compare.haplotypes[!is.na(haplotype), unitig_name]
d <- tibble(value=unitigs.set, `haplotagged unitigs`=haplotagged.unitigs, `haplo-clustered unitigs`=haploclust.unitigs)
venn.diagram <- ggplot(d)+geom_venn(aes(A=`haplotagged unitigs`,B=`haplo-clustered unitigs`))+coord_fixed()+theme_void()
venn.diagram

# haploclust accuracy
compare.haplotypes[, `:=`(same.haplo=0, diff.haplo=0)]#, diff.chrom=0)]
#compare.haplotypes[ref_name!=chrom_clust, diff.chrom:=1]
compare.haplotypes[ref_name==chrom_clust & haplotype==unitig_gr_haplo, `:=`(same.haplo=1, diff.haplo=0)]
compare.haplotypes[ref_name==chrom_clust & haplotype!=unitig_gr_haplo, `:=`(same.haplo=0, diff.haplo=1)]
#compare.haplotypes[haplotype==unitig_gr_haplo, `:=`(same.haplo=1, diff.haplo=0)]
#compare.haplotypes[haplotype!=unitig_gr_haplo, `:=`(same.haplo=0, diff.haplo=1)]

compare.haplotypes[, `:=`(num.haplo.match=sum(same.haplo), num.haplo.mismatch=sum(diff.haplo)), by=ref_name]
# For chroms that have haplotype switch
haplo.set = c('H1','H2')
compare.haplotypes[num.haplo.match<num.haplo.mismatch & unitig_gr_haplo=='H1', unitig_gr_haplo:='2']
compare.haplotypes[num.haplo.match<num.haplo.mismatch & unitig_gr_haplo=='H2', unitig_gr_haplo:='1']
compare.haplotypes[unitig_gr_haplo %in% c('1','2'), unitig_gr_haplo:=paste0('H', unitig_gr_haplo)]

compare.haplotypes[, prediction_status:='homozygous']
compare.haplotypes[ref_name!=chrom_clust, prediction_status:='false chrom cluster']
compare.haplotypes[ref_name==chrom_clust & unitig_gr_haplo!=haplotype, prediction_status:='false haplo cluster']
compare.haplotypes[ref_name==chrom_clust & is.na(unitig_gr_haplo) & !is.na(haplotype), prediction_status:='false haplo cluster']
compare.haplotypes[ref_name==chrom_clust & unitig_gr_haplo==haplotype, prediction_status:='true chrom/haplo cluster']
compare.haplotypes[ref_name==chrom_clust & !is.na(unitig_gr_haplo) & is.na(haplotype), prediction_status:='not clustered by haploclust']
# the label is empty string iff the bot predicted and gr haplo are NA

accuracy <- compare.haplotypes[!is.na(unitig_gr_haplo), .(ref_name, same.haplo, diff.haplo, num.haplo.match, num.haplo.mismatch)]
accuracy[,`:=`(true_prediction=max(num.haplo.match,num.haplo.mismatch), false_prediction=min(num.haplo.match,num.haplo.mismatch)), by=ref_name]
accuracy <- accuracy[, head(.SD, 1), by=ref_name]
accuracy[, acc:=true_prediction/(true_prediction+false_prediction)]

true.total <- accuracy[,sum(true_prediction)]
false.total <- accuracy[,sum(false_prediction)]
message('whole-genome haplotype clustering accuracy = ', true.total/(true.total+false.total))

# plotting barplot prediction type per chromosome
# counting prediction status by chromosome
prediction.labels <- c('homozygous','false chrom cluster','false haplo cluster','not clustered by haploclust', 'true chrom/haplo cluster')
prediction.colors <- c('homozygous'='white', 'false chrom cluster'='red', 'false haplo cluster'='yellow', 'not clustered by haploclust'='gray', 'true chrom/haplo cluster'='blue')

counts <- compare.haplotypes %>% group_by(ref_name) %>% count(prediction_status) %>% mutate(percent=n/sum(n))
counts <- as.data.table(counts)
ord.chr <- factor(counts[, ref_name], levels=valid.chroms)
ord.pred.status <- factor(counts[, prediction_status], levels=prediction.labels)
ggplot(counts, aes(x=ord.chr, y=percent, fill=ord.pred.status))+
  geom_bar(stat='identity', width = 0.75, color='black')+scale_y_continuous(labels=scales::percent)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  xlab('chromosome')+ylab('relative frequencies')+ylim(0,1)+
  scale_fill_manual('prediction status', values=prediction.colors)

# pie chart
total.counts <- counts[, .(percent=sum(percent)), by=prediction_status]
total.counts[, percent:=percent/sum(percent)]

# Compute the position of labels
total.counts <- total.counts %>%
  arrange(desc(prediction_status)) %>%
  mutate(ypos = cumsum(percent)- 0.5*percent )
#pi.chart <- 
ggplot(total.counts, aes(x='', y=percent, fill=prediction_status))+
  geom_bar(stat='identity', width = 1, color='black')+scale_y_continuous(labels=scales::percent)+coord_polar("y", start=0)+
  scale_fill_manual('prediction status', values=prediction.colors)+
  geom_text(aes(y = ypos, label = percent(percent)), color = "black", size=3)+
  theme_void()
