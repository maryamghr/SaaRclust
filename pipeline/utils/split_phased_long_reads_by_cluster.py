def print_dict_head(dic, num):
	for i in range(num):
		key = list(dic.keys())[i]
		print(key, ':', dic[key])


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
				name, clust = sp[0], sp[1]

				# remove the ref alignment info from read name
				name_sp = name.split('/ccs')
				name = name_sp[0]+'/ccs'

				read_to_clust[name] = clust

	return read_to_clust


def split_phase_files_by_clust(long_reads_phase_files, long_reads_clust_files, output_files):

	cluster_to_output_file = {}
	read_to_clust = get_cluster(long_reads_clust_files)

	print('head read_to_clust:')
	print_dict_head(read_to_clust,5)
	
	for out in output_files:
		base_name = out.split('/')[-1]
		cluster = base_name.split('_')[0]
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
				name = name_sp[0]+'/ccs'

				if haplo=='0':
					haplo='H1'
				elif haplo=='1':
					haplo='H2'
				else:
					haplo='none'
				
				if name not in read_to_clust:
					# the read is not clustered
					continue

				clust = read_to_clust[name]

				if clust not in cluster_to_output_file:
					# the cluster is the garbage cluster
					continue

				out_file = cluster_to_output_file[clust]
				out_file.write(name + '\t' + haplo + '\n')

	for clust in cluster_to_output_file:
		out_file = cluster_to_output_file[clust]
		out_file.close()
