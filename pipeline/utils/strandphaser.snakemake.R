log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source('utils/bubble_phasing_lts.R')
source('utils/process_mummer_map.R')

sample=snakemake@wildcards[['sample']]
print(paste('sample=', sample))

clusters <- strsplit(snakemake@wildcards[["clust_pair"]], "_")[[1]]
clusters <- c(min(clusters), max(clusters))

print('got clusters')

wc.cell.clust <- fread(snakemake@input[["wc_cell_clust"]])
ss.clust <- fread(snakemake@input[["ss_clust"]], header=F)

print('got ss clust')
print(paste('map:', snakemake@input[["map"]]))
print(paste('bubbles:', snakemake@input[["bubbles"]]))
#map <- output_mummer_map_table(snakemake@input[["map"]], snakemake@input[["bubbles"]])
#map <- map[!is.na(SSstart)]

#print('got mummer map table')

# outputting map tables
#fwrite(map, file=snakemake@output[["valid_maps"]], row.names=F, sep='\t')

map <- fread(snakemake@input[["map"]])
map <- map[!is.na(bubbleAllele)]
map.sp <- output_bubble_allele_coverage_matrix(clusters[1], clusters[2], wc.cell.clust, ss.clust, map)

print('splitted map')

## Get selected library names
select.libs <- unique(wc.cell.clust$lib)

print('select.libs')
print(select.libs)

strandphaser(map.sp[[clusters[1]]], map.sp[[clusters[2]]], snakemake@wildcards[["clust_pair"]], select.libs, snakemake@output[["phased_strand_states"]])

print ('done')
