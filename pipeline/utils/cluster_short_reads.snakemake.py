from cluster_short_reads import *
import time

clust_cov_files = snakemake.input['clust_cov']
log_file = snakemake.log[0]
colnames = ['clust.forward', 'name']
outfile = snakemake.output[0]

print('clust_cov_files:')
print(clust_cov_files)
print('log file =', log_file)

print('start clustering...')

with open(log_file, 'w') as log:
	print('clustering short reads', file=log)
	start_time = time.time()
	short_read_to_clust_cov = cluster_short_reads_reads(clust_cov_files)
	output_short_reads_clust(short_read_to_clust_cov, outfile, colnames)
	print('elapsed time =', time.time(), file=log)

