import gzip
import time


def cluster_ss_reads(ss_clust_cov_files):
	print('cluster_ss_reads2')
	ss_clust_to_cov = {}

	for cov_file in ss_clust_cov_files:
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

				if (ss, clust) not in ss_clust_to_cov:
					ss_clust_to_cov[(ss,clust)] = cov
				else:
					ss_clust_to_cov[(ss,clust)] += cov
		print('elapsed time =', time.time()-start_time, 's')

	# assiging each ss read to the cluster with maximum coverage

	ss_to_clust_cov = {}

	print('assigning SS reads to clusters...')
	for (ss, clust) in ss_clust_to_cov:
		cov = ss_clust_to_cov[(ss, clust)]

		if ss not in ss_to_clust_cov:
			ss_to_clust_cov[ss] = (clust, cov)
		else:
			curr_ss_cov = ss_to_clust_cov[ss][1]
			if curr_ss_cov < cov:
				ss_to_clust_cov[ss] = (clust, cov)


	return ss_to_clust_cov


def output_ss_clust(ss_to_clust_cov, out_file, colnames):
	with open(out_file, 'w') as out:
		print(colnames[0], '\t', colnames[1], file=out)
		for ss in ss_to_clust_cov:
			print(ss_to_clust_cov[ss][0], '\t', ss, file=out)



