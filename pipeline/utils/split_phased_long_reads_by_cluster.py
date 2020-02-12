

def get_cluster(long_reads_clust_files):
	read_to_clust = {}
	for clust_file in long_reads_clust_files:
		with open(clust_file) as f:
			# skip header
			next(f)

			for line in f:
				if line=="":
					continue

				sp = line.split()
				read_to_clust[sp[0]] = sp[1]

	return read_to_clust


def split_phase_files_by_clust(long_reads_phase_files, long_reads_clust_files, output_files):

	cluster_to_output_file = {}
	read_to_clust = get_cluster(long_reads_clust_files)
	
	for out in output_files:
		cluster = out.split('_')[0]
		# remove 'cluster' from the cluster name
		cluster = cluster.split('cluster')[1]
		out_file = open(out, 'w')
		cluster_to_output_file[cluster] = out_file
		

	for phase_file in long_reads_phase_files:
		print('reading', phase_file)
		with open(phase_file) as f:
			# skip header
			next(f)

			for line in f:
				if line=="":
					continue

				sp = line.split()

				name, haplo = sp[0], sp[-1]
				# remove the ref alignment info from read name
				name_sp = name.split('/ccs')
				name = name_sp[0][1:]+'/ccs'
				print('name =', name)

				if name not in read_to_clust:
					# the read is not clustered
					continue

				clust = read_to_clust[name]
				out_file = cluster_to_output_file[clust]
				out_file.write(name + '\t' + haplo + '\n')

	for clust in cluster_to_output_file:
		out_file = cluster_to_output_file[clust]
		out_file.close()
