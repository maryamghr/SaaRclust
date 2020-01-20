import gzip
import time

def map_cluster_to_pair_and_chrom_flag(clust_partners_file):
	cluster_to_chrom_flag = {}
	cluster_pair = {}
	with open(clust_partners_file) as f:
		# skip the header line
		next(f)

		for line in f:
			if line=="":
				break
			sp = line.split()
			cluster_to_chrom_flag[sp[1]]=sp[0]
			cluster_pair[sp[1]]=sp[2]
			
	return cluster_to_chrom_flag, cluster_pair
	
def cluster_ss_reads(minimap_files, cluster_pair, out_file, colnames):
	ss_clust_count = {}
	ss_chrom_flag = {}
	ss_to_clust = {}
	for minimap in minimap_files:
		start_time = time.time()
		print('processing', minimap)
		with gzip.open(minimap) as f:
			# skip the header line
			next(f)
			for line in f:
				line = line.decode()
				if line=="":
					break

				sp = line.split()

				ss_name, ss_flag, ss_chrom, strand, clust = sp[1], sp[3], sp[4], sp[9], sp[19]
				#print('ss_name, ss_flag, ss_chrom, strand, clust:')
				#print(ss_name, ss_flag, ss_chrom, strand, clust)

				if clust not in cluster_pair:
					# the cluster does not have a pair (garbage cluster)
					#print(clust, 'does not exist in cluster pairs')
					continue

				if ss_name not in ss_chrom_flag:
					ss_chrom_flag[ss_name] = ss_chrom + '_' + ss_flag

				if strand == '-':
					clust = cluster_pair[clust]

				
				# adding the count of cluster for the ss name
				if ss_name not in ss_clust_count:
					ss_clust_count[ss_name] = {clust:1}

				elif clust not in ss_clust_count[ss_name]:
					ss_clust_count[ss_name][clust] = 1

				else:
					ss_clust_count[ss_name][clust] += 1
		print('elapsed time =', time.time()-start_time, 's')

	print('calling the clusters with max count for every ss read...')

	with open(out_file, 'w') as out:
		print(colnames[0], '\t', colnames[1], file=out)
		for ss in ss_clust_count:
			max_count, clust_with_max_count = 0, 'None'
			for clust in ss_clust_count[ss]:
				count = ss_clust_count[ss][clust]
				if count > max_count:
					max_count, clust_with_max_count = count, clust

			ss_to_clust[ss] = clust_with_max_count
			print(clust_with_max_count, '\t', ss, file=out)

	return (ss_to_clust, ss_chrom_flag)


def evaluate_clustering_ss_reads(ss_to_clust, ss_chrom_flag, cluster_to_chrom_flag, out_file):
	print('evaluating clutering ss reads')
	with open(out_file, 'w') as out:
		print(len(ss_to_clust), 'ss reads are clustered.', file=out)
		
		true, false = 0, 0
		for ss in ss_to_clust:
			clust = ss_to_clust[ss]
			if cluster_to_chrom_flag[clust] == ss_chrom_flag[ss]:
				true += 1
			else:
				false += 1

		print('#correctly clustered ss reads =', true, file=out)
		print('#wrongly clustered ss reads =', false, file=out)




#start_time = time.time()
#minimap_file = pd.read_csv('aligns_k15_w1_f0.1_z500/HG00514_chunk000_withclust.maf.gz', usecols=['SSreadNames','SSflag', 'SSchrom', 'strand', 'ML_cluster_idx'], compression='gzip', sep='\t')
#print('elapsed time =', time.time()-start_time, 's')


