library(data.table)
library(seqinr)
library(Rsamtools)
library(ggplot2)

#setwd('/home/maryam/research/haploclust/haploclust/pipeline/')
bam.file <- "../../HG00733/hifiasm/ref_aln/asm.r_utg.haplotagged.bam"
predicted.haplo.files <- list.files('../../HG00733/phased_unitigs/', full.names = T)
predicted.chrom.file <- '../../HG00733/SaaRclust/Clusters/MLclust.data'
valid.chroms <- paste0('chr', c(1:22,'X'))

unitigs.bam <- scanBam(BamFile(bam.file), param=ScanBamParam(what=c("rname", "qname"), tag="HP", flag=scanBamFlag(isSupplementaryAlignment=FALSE)))
original.haplotypes <- data.table(ref_name=unitigs.bam[[1]]$rname, unitig_name=unitigs.bam[[1]]$qname, unitig_gr_haplo=as.character(unlist(unitigs.bam[[1]]$tag)))
original.haplotypes[!is.na(unitig_gr_haplo), unitig_gr_haplo:=paste0('H',unitig_gr_haplo)]
predicted.haplotypes <- lapply(predicted.haplo.files, fread)
predicted.haplotypes <- Reduce(rbind, predicted.haplotypes)
names(predicted.haplotypes)[1] <- 'unitig_name'
predicted.chroms <- fread(predicted.chrom.file)
names(predicted.chroms)[1] <- 'unitig_name'

compare.haplotypes <- merge(original.haplotypes, predicted.haplotypes, by='unitig_name')
compare.haplotypes <- merge(compare.haplotypes, predicted.chroms, by='unitig_name', all.x=T)
#compare.haplotypes[, opposite_haplotype:=""]
#compare.haplotypes[haplotype=='H2', opposite_haplotype:='H1']
#compare.haplotypes[haplotype=='H1', opposite_haplotype:='H2']

compare.haplotypes[, `:=`(same.haplo=0, diff.haplo=0)]
compare.haplotypes[haplotype==unitig_gr_haplo, `:=`(same.haplo=1, diff.haplo=0)]
compare.haplotypes[haplotype!=unitig_gr_haplo, `:=`(same.haplo=0, diff.haplo=1)]

accuracy <- compare.haplotypes[, .(ref_name, same.haplo, diff.haplo)]
accuracy[, `:=`(num.haplo.match=sum(same.haplo), num.haplo.mismatch=sum(diff.haplo)), by=ref_name]
# For chroms that have haplotype switch
accuracy[,`:=`(true_prediction=max(num.haplo.match,num.haplo.mismatch), false_prediction=min(num.haplo.match,num.haplo.mismatch)), by=ref_name]
accuracy <- accuracy[ref_name %in% valid.chroms, head(.SD, 1), by=ref_name]
accuracy[, acc:=true_prediction/(true_prediction+false_prediction)]

true.total <- accuracy[,sum(true_prediction)]
false.total <- accuracy[,sum(false_prediction)]
message('whole-genome haplotype clustering accuracy = ', true.total/(true.total+false.total))

# plotting the accuracy
acc.plt <- ggplot(accuracy, aes(x=ref_name,y=acc,fill='red'))+geom_bar(stat='identity', width = 0.75)