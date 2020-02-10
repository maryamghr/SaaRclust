import gzip
import time


def cluster_short_reads_reads(short_reads_clust_cov_files):
	print('cluster_short_reads_reads')
	short_reads_clust_to_cov = {}

	for cov_file in short_reads_clust_cov_files:
		print('processing file', cov_file)
		start_time = time.time()
		with open(cov_file) as f:
			# skip the header line
			next(f)

			for line in f:
				if line == "":
					break

				sp = line.split()
				ss, clust, cov = sp[0], sp[1], int(sp[2])

				if (ss, clust) not in short_reads_clust_to_cov:
					short_reads_clust_to_cov[(ss,clust)] = cov
				else:
					short_reads_clust_to_cov[(ss,clust)] += cov
		print('elapsed time =', time.time()-start_time, 's')

	# assiging each ss read to the cluster with maximum coverage

	short_reads_to_clust_cov = {}

	print('assigning short reads to clusters...')
	for (ss, clust) in short_reads_clust_to_cov:
		cov = short_reads_clust_to_cov[(ss, clust)]

		if ss not in short_reads_to_clust_cov:
			short_reads_to_clust_cov[ss] = (clust, cov)
		else:
			curr_short_reads_cov = short_reads_to_clust_cov[ss][1]
			if curr_short_reads_cov < cov:
				short_reads_to_clust_cov[ss] = (clust, cov)


	return short_reads_to_clust_cov


def output_short_reads_clust(short_reads_to_clust_cov, out_file, colnames):
	with open(out_file, 'w') as out:
		print(colnames[0], '\t', colnames[1], file=out)
		for ss in short_reads_to_clust_cov:
			print(short_reads_to_clust_cov[ss][0], '\t', ss, file=out)



