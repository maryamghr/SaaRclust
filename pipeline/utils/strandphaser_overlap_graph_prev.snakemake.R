log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

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

bubble.cov.files <- unlist(snakemake@input[["bubble_lib_cov_matrix"]])

print('bubble.cov.files =')
print(bubble.cov.files)

clust.pairs <- read.table(snakemake@input[["clust_to_chrom"]], header = TRUE)

## Keep only unique pairs
forw <- as.numeric(gsub("[^\\d]+", "", clust.pairs$clust.forward, perl=TRUE))
backw <- as.numeric(gsub("[^\\d]+", "", clust.pairs$clust.backward, perl=TRUE))
pairs <- cbind(pmin(forw, backw), pmax(forw, backw))
clust.pairs <- clust.pairs[!duplicated(pairs),]

print('clust.pairs =')
print(clust.pairs)

output_phased_strand_states(bubble.cov.files, clust.pairs, select.libs, snakemake@output[["phased_strand_states"]], snakemake@output[["phased_bubbles"]])
