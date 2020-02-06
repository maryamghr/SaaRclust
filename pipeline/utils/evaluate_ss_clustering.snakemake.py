from evaluate_ss_clustering import *

print('getting clust_to_chrom_dir')
clust_to_chrom_dir = get_clust_to_chrom_dir(snakemake.input['clust_to_chrom_dir'])
print('getting ss_to_clust')
ss_to_clust = get_ss_clusters(snakemake.input['ss_clust'])
print('getting ss_to_chrom_dir')
ss_to_chrom_dir = get_ss_chrom_dir(snakemake.input['ss_bam_list'])
print('evaluating ss clustering')
evaluate_ss_clustering(ss_to_chrom_dir, ss_to_clust, clust_to_chrom_dir, snakemake.output[0])
