log <- file(snakemake@log[[1]], open='wt')
#sink(file=log, type='message')
#sink(file=log, type='output')

source('utils/bubble_phasing_lts.R')
sample=snakemake@wildcards[['sample']]

print(paste('sample=', sample))

WC.regions <- read.table(snakemake@input[["wc_cell_clust"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

print('WC.regions')
print(WC.regions)
## Get selected library names
select.libs <- unique(WC.regions$lib)

print('select.libs')
print(select.libs)


strandphaser(snakemake@input[["clust1"]], snakemake@input[["clust2"]], snakemake@wildcards[["clust_pair"]], select.libs, snakemake@output[["phased_strand_states"]]) #, snakemake@output[["phased_bubbles"]])


