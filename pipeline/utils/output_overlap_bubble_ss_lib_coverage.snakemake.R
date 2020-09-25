log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(data.table)
library(reshape2)
library(seqinr)
source('utils/process_mummer_map.R')

print('snakemake@input[["valid_maps"]]:')
print(snakemake@input[["valid_maps"]])
print('snakemake@input[["ss_clust"]]')
print(snakemake@input[["ss_clust"]])
print('snakemake@output:')
print(snakemake@output)
print("snakemake@output[[\"clust1\"]]")
print(snakemake@output[["clust1"]])
print("snakemake@output[[\"clust2\"]]")
print(snakemake@output[["clust2"]])


wc.cell.clust <- fread(snakemake@input[["wc_cell_clust"]])
ss.clust <- fread(snakemake@input[["ss_clust"]], header=F)

map <- fread(snakemake@input[["valid_maps"]])
map.sp <- output_bubble_allele_coverage_matrix <- function(snakemake@wildcards[["clust_pair"]], wc.cell.clust, ss.clust, map)

fwrite(map.sp[[clusters[1]]], file=snakemake@output[["clust1"]], sep="\t")
fwrite(map.sp[[clusters[2]]], file=snakemake@output[["clust2"]], sep="\t")


# run in R console:
# wc.cell.clust <- fread('aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/wc_cells_clusters.data')
# ss.clust <- fread('aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/split/ss_reads/clusterV25_V32.data', header=F)
# map <- fread('aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/split/hifiasm/bubble_ss_map/valid_V25_V32_r_utg_maximal_uniqe_exact_match.data')
# clusters <- strsplit("V25_V32", "_")[[1]]
