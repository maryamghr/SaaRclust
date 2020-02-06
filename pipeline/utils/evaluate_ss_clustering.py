from __future__ import division
import pysam
import time

def get_clust_to_chrom_dir(clust_to_chrom_dir_file):
	clust_to_chrom_dir = {}
	with open(clust_to_chrom_dir_file) as f:
		for line in f:
			if not line.startswith('chr'):
				# it is a hedear line
				continue

			chrom_dir, clust_forward, clust_backward = line.split()

			if clust_forward not in clust_to_chrom_dir:
				# remove the direction from chrom name
				clust_to_chrom_dir[clust_forward] = chrom_dir

	return clust_to_chrom_dir



def get_ss_clusters(ss_clust_file):
	ss_to_clust = {}
	with open(ss_clust_file) as f:
		#skip the header line
		next(f)

		for line in f:
			if line=="":
				continue

			ss_clust, ss_name = line.split()

			ss_to_clust[ss_name] = ss_clust

	return ss_to_clust


def get_ss_chrom_dir(ss_bam_list):
	ss_to_chrom_dir = {}
	valid_chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX']
	
	for ss_bam in ss_bam_list:
		start_time = time.time()
		print('processing', ss_bam, '...')
		bamfile = pysam.AlignmentFile(ss_bam, 'rb')

		for read in bamfile.fetch():
			chrom = read.reference_name

			if chrom not in valid_chroms:
				continue

			direction = '16' if read.is_reverse else '0'
			ss_name = read.query_name
			
			ss_to_chrom_dir[ss_name] = chrom + '_' + direction

		print('elapsed time =', time.time()-start_time, 's')

	return ss_to_chrom_dir
	
	
def evaluate_ss_clustering(ss_to_chrom_dir, ss_to_clust, clust_to_chrom_dir, out_file):
	with open(out_file, 'w') as out:
		print('total number of ss reads with valid chromsome/direction =', len(ss_to_chrom_dir), file=out)
		print('total number of clustered ss reads =', len(ss_to_clust), file=out)
		
		true, false = 0, 0
		for ss_name in ss_to_clust:
			ss_clust = ss_to_clust[ss_name]
			
			if ss_clust not in clust_to_chrom_dir:
				continue
				
			if clust_to_chrom_dir[ss_clust] == ss_to_chrom_dir[ss_name]:
				true += 1
			else:
				false+=1
				
		print('total number of clustered ss reads with non-garbage(non-repeat) clusters =', true+false, file=out)
		print('number of true clustered ss reads =', true, ',' , true*100/(true+false), '% of clustered reads', file=out)
		print('number of false clustered ss reads =', false, ',' , false*100/(true+false), '% of clustered reads', file=out)
		



