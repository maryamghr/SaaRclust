log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)

d = fread(snakemake@input[[1]])
trim = snakemake@params[["trim"]]
br = snakemake@params[["breaks"]]

tr = quantile(d$V1, probs = c(trim, 1-trim))
pdf(snakemake@output[[1]])
hist(d[V1>tr[1] & V1<tr[2], V1], breaks=br, xlab="kmer coverage", ylab="kmer count", main="contigs kmer-coverage histogram")
dev.off()
