from cluster_ss_reads import *

clust_partners_file = snakemake.input['clust_partners']
minimap_files = snakemake.input['minimap']
out_file = snakemake.output[0]
log_file = snakemake.log[0]
colnames = ['clust.forward', 'name']

cluster_to_chrom_flag, cluster_pair = map_cluster_to_pair_and_chrom_flag(clust_partners_file)
print(cluster_to_chrom_flag)
print(cluster_pair)
ss_to_clust, ss_chrom_flag = cluster_ss_reads(minimap_files, cluster_pair, out_file, colnames)
evaluate_clustering_ss_reads(ss_to_clust, ss_chrom_flag, cluster_to_chrom_flag, log_file)
