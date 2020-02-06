from cluster_ss_reads2 import *
import time

ss_clust_cov_files = snakemake.input['ss_clust_cov']
log_file = snakemake.log[0]
colnames = ['clust.forward', 'name']
outfile = snakemake.output[0]

print('ss_clust_cov_files:')
print(ss_clust_cov_files)
print('log file =', log_file)

print('start clustering...')

with open(log_file, 'w') as log:
	print('clustering ss reads', file=log)
	start_time = time.time()
	ss_to_clust_cov = cluster_ss_reads(ss_clust_cov_files)
	output_ss_clust(ss_to_clust_cov, outfile, colnames)
	print('elapsed time =', time.time(), file=log)

