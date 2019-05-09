log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

print(snakemake@output[["acc_plot"]])
print(snakemake@output[["acc_table"]])
print(snakemake@params[['minLib']])

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), 'utils/R-packages/'))

suppressPackageStartupMessages(library(SaaRclust))
library(data.table)


source('utils/SaaRclust_evaluation_plots.R')

inputfolder <- snakemake@input[['inputfolder']]
thresholds <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
minLib <- snakemake@params[['minLib']]

plt <- ClustersAccuracyPerChrPerDir(inputfolder, thresholds, minLib)

ggsave(snakemake@output[["acc_plot"]], plt$acc.plot)
fwrite(plt$plot.table, snakemake@output[["acc_table"]], sep="\t")
