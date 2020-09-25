log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)

phase <- fread(snakemake@input[["phase"]])
soft.clust <- fread(snakemake@input[["clust"]], header=F, col.names='long_read_name')

clust.phase <- merge(soft.clust, phase, by='long_read_name')
clust.phase[haplotype=='None', haplotype:="none"]
clust.phase[haplotype=='0', haplotype:="H1"]
clust.phase[haplotype=='1', haplotype:="H2"]

fwrite(clust.phase, snakemake@output[[1]], sep="\t", row.names=F, col.names=F, quote=F)
