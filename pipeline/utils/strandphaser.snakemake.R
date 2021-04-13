log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source('utils/bubble_phasing_lts.R')
source('utils/process_mummer_map.R')

sample=snakemake@wildcards[['sample']]
print(paste('sample=', sample))

clust <- snakemake@wildcards[["clust"]]
clust.pairs <- fread(snakemake@input[["clust_pairs"]])
clust.pairs <- clust.pairs[chrom_clust==clust]
clusters <- c(clust.pairs$first_clust, clust.pairs$second_clust)
clusters <- as.character(clusters)
wc.cells.clust <- fread(snakemake@input[["wc_cell_clust"]])

print('got clusters')

wc.cell.clust <- fread(snakemake@input[["wc_cell_clust"]])
ss.clust <- fread(snakemake@input[["ss_clust"]], header=F)

print(ss.clust)
print(paste('map:', snakemake@input[["map"]]))
print(paste('bubbles:', snakemake@input[["bubbles"]]))


map <- fread(snakemake@input[["map"]])
map <- map[bubbleAllele!="None"]
map.sp <- output_bubble_allele_coverage_matrix(clusters, wc.cell.clust, ss.clust, map)

print('splitted map')

## Get selected library names
select.libs <- wc.cells.clust[clust.forward %in% clusters, unique(lib)]

print('select.libs')
print(select.libs)

strandphaser(map.sp[[clusters[1]]], map.sp[[clusters[2]]], snakemake@wildcards[["clust"]], select.libs, snakemake@output[["phased_strand_states"]])

print ('done')
