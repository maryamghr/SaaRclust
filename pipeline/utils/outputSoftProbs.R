sink(snakemake@log[[1]])
print(paste("input=", snakemake@input[[1]]))
print(paste("output=", snakemake@output[[1]]))
library(data.table)
load(snakemake@input[[1]])
dt <- data.table(PBname=rownames(soft.clust.obj$soft.pVal),
                 PBchrom=soft.clust.obj$PBchrom,
		 PBflag=soft.clust.obj$PBflag,
                 soft.clust.obj$soft.pVal)

write.table(dt, file=snakemake@output[[1]], sep="\t", quote = F, row.names = F)
