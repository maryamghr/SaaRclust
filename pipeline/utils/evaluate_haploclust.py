from __future__ import division
from bubble_long_read_alignment import *
import pysam
import time
import pdb
from argparse import ArgumentParser

global valid_chroms
valid_chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX']

def print_dict_head(dic, num):
	for i in range(num):
		key = list(dic.keys())[i]
		print(key, ':', dic[key])




def get_bubbles(bubble_haplotagged_bam_file, with_km=True):

	'''
	Reads bubbles happlotagged bam files and returns the list of bubbles for each valid chromosome

	Parameters (list(str)):
		bubble_haplotagged_bam_file: A list of paths of bubbles haplotagged bam files

	Returns (dict(int): {bubble_id -> bubble}):
		A dictionary that maps each bubble_id to its bubble object
		
	'''

	start_time = time.time()
	print('getting bubbles...')

	bubbles = {}
	
	samfile = pysam.AlignmentFile(bubble_haplotagged_bam_file, 'rb')

	for read in samfile.fetch(until_eof=True): # until_eof=True allows to read also unmapped reads in sam file
		bubble_name = read.query_name
		bubble_name_sp = bubble_name.split('_')
		bubble_id, bubble_allele_id = int(bubble_name_sp[1]), int(bubble_name_sp[3])-1

		if bubble_id in bubbles:
			bubble = bubbles[bubble_id]

		else:
			bubble = Bubble(bubble_id)
			bubbles[bubble_id] = bubble

		bubble_allele = BubbleAllele(bubble_allele_id, bubble)
		bubble.add_allele(bubble_allele)

		if with_km:
			bubble_allele.km = float(bubble_name_sp[7])
		
		if read.is_unmapped:
			bubble.actual_type = 'unmapped'
			continue


		chrom = read.reference_name

		if not chrom in valid_chroms:
			# the chromosome name is not valid
			bubble.actual_type = 'invalid_chrom'
			continue		

		bubble.actual_chrom = chrom

		if not read.has_tag("HP"):
			# the read is not haplotagged
			bubble.actual_type = 'untagged'
			continue

		bubble.actual_type = 'tagged' + str(read.get_tag("HP")-1)
		bubble_allele.actual_haplo = read.get_tag("HP")-1

	print('num bubbles =', len(bubbles))
	print('after the for loop', bubbles[bubble_id])

	print('elapsed time =', time.time()-start_time)		
		
	return bubbles


def get_clust_to_chrom(clust_to_chrom_file):

	'''
	'''

	start_time = time.time()
	clust_to_chrom = {}
	with open(clust_to_chrom_file) as f:
		for line in f:
			if not line.startswith('chr'):
				# it is a hedear line
				continue

			chrom, clust_forward, clust_backward = line.split()

			if clust_forward not in clust_to_chrom:
				# remove the direction from chrom name
				clust_to_chrom[clust_forward] = chrom.split('_')[0]

	print('elapsed time =', time.time()-start_time)	

	return clust_to_chrom


# from now on, change chrom_bubbles to bubbles
# TODO: change it later in the get_bubbles file...

def add_bubble_clust(bubble_clust_file, bubbles):

	'''
	Adds SaaRclust (chrom+dir) clusters to bubble objects

	Parameters (str):
		bubble_clust_file (str): Path of bubbles clusters file

		bubbles (dict(int): {bubble_id -> bubble}):
		A dictionary that maps each bubble_id to its bubble object
		
	'''

	start_time = time.time()
	print('adding SaaRclust chrom clusters to bubbles')

	bubble_id_to_clust = {}
	with open(bubble_clust_file) as f:
		for line in f:
			sp = line.split()

			if not sp[0].startswith('V'):
				# it is a header line
				continue

			bubble_id, bubble_clust = int(sp[1]), sp[0]

			assert (bubble_id in bubbles), 'bubble ' + str(bubble_id) + ' is not present in the bubbles'

			#if bubble_id not in bubbles:
				# bubble is not present in the bubbles (not mapped to a valid chromosome)
			#	continue

			bubble = bubbles[bubble_id]
			bubble.clust = bubble_clust

	print('elapsed time =', time.time()-start_time)	


def add_bubble_allele_pred_haplo(bubble_phase_file, bubbles):

	'''
	Adds Haploclust predicted haploypes to bubble allele objects

	Parameters (str):
		bubble_phase_file (str): Path of bubbles phase file

		bubbles (dict(int): {bubble_id -> bubble}):
		A dictionary that maps each bubble_id to its bubble object
		
	'''

	start_time = time.time()
	print('adding haplotype clusers to bubble alleles')

	with open(bubble_phase_file) as f:
		# skip the header line
		next(f)

		for line in f:
			sp = line.split()
			bubble_id, al0_haplo = int(sp[0]), int(sp[1])

			assert (bubble_id in bubbles), 'bubble ' + str(bubble_id) + ' is not present in the bubbles'

			#if bubble_id not in bubbles:
				# bubble is not present in the bubbles (not mapped to a valid chromosome)
			#	continue

			bubble = bubbles[bubble_id]

			bubble.allele0.pred_haplo = al0_haplo
			bubble.allele1.pred_haplo = 1-al0_haplo

	print('elapsed time =', time.time()-start_time)


def get_long_reads(long_reads_haplotagged_bam_files):

	'''
	'''

	start_time = time.time()
	print('getting long reads')

	long_reads = {}	
	
	for alignmentfile in long_reads_haplotagged_bam_files:
		print('processing', alignmentfile, '...')
		samfile = pysam.AlignmentFile(alignmentfile, 'rb')
		for read in samfile.fetch(until_eof=True):
			read_name = read.query_name
		
			
			long_read = LongRead(read_name)
			long_reads[read_name] = long_read

			if read.is_unmapped:
				long_read.actual_type = 'unmapped'
				continue

			chrom = read.reference_name
			
			if not chrom in valid_chroms:
				# the chromosome name is not valid
				long_read.actual_type = 'invalid_chrom'
				continue		
			
			long_read.actual_chrom = chrom

			if not read.has_tag("HP"):
				# the read is not haplotagged
				long_read.actual_type = 'untagged'
				continue

			haplotype = read.get_tag("HP")-1
			long_read.actual_haplo = haplotype
			long_read.actual_type = 'tagged' + str(haplotype)
			
	return long_reads 


def add_long_reads_pred_haplotype(long_reads_phase_file_list, long_reads):

	'''
	Adds Haploclust predicted haploypes to bubble allele objects

	Parameters (str):
		long_reads_phase_file_list (str): A list of paths of bubbles haplotagged bam files

		long_reads (dict(str): {long_read name -> long_read}):
		A dictionary that maps each read name to its LongRead object
		
	'''

	start_time = time.time()
	print('adding pred haplo to long reads')

	for long_reads_phase_file in long_reads_phase_file_list:
		with open(long_reads_phase_file) as f:
			# skip the header line
			next(f)
			for line in f:
				sp = line.split()
				read_name, haplo = sp[0], sp[-1]
				read_name_sp = read_name.split('/ccs')
				read_name = read_name_sp[0]+'/ccs'

				assert (read_name in long_reads), 'read ' + read_name + ' should be present in long_reads'

				long_read = long_reads[read_name]

				if haplo!="?":
					long_read.pred_haplo = int(haplo)

				long_read.haplo0_edit_dist, long_read.haplo1_edit_dist = (int(sp[3]), int(sp[4]))


	print('elapsed time =', time.time()-start_time)


def set_alignments(kmers_files_list, bubbles, phased_long_reads):

	'''
	'''

	start_time = time.time()
	print('getting alignments')

	alignments = {}

	for kmers_files in kmers_files_list:
		print('reading kmers from file', kmers_files)
		with open(kmers_files) as f:
			# skip the header line
			next(f)
			for line in f:
				sp = line.split()
				if sp[0]=='none':
					continue # the long read is not mapped to any bubble

				bubble_id, bubble_al, read_name, edit_dist = \
				int(sp[0]), int(sp[1]), sp[2], int(sp[5])

				read_name_sp = read_name.split('/ccs')
				read_name = read_name_sp[0]+'/ccs'

				assert (bubble_id in bubbles), 'bubble ' + str(bubble_id) + ' is not present in the bubbles'
				assert (read_name in long_reads), 'long read ' + read_name + ' is not present in long reads'

				#if bubble.pred_haplo == None or read_name not in phased_long_reads:
				#	# either the bubble or the read_name is not phased
				#	continue

				bubble = bubbles[bubble_id]
				long_read = phased_long_reads[read_name]

				bubble_allele = bubble.allele0 if bubble_al==0 else bubble.allele1

				assert (bubble_allele != None), 'bubble ' + str(bubble_id) + ' allele ' + str(bubble_al) + ' is None'
				
				aln = Alignment(long_read, bubble_allele, edit_dist)

	print('elapsed time =', time.time()-start_time)	


def evaluate_bubble_clustering(bubbles, clust_to_chrom, output_file):

	'''
	'''

	start_time = time.time()
	print('evaluating bubbles clustering')

	num_bubbles=len(bubbles)
	num_chr_clustered_bubbles=0
	num_garbage_chr_clustered_bubbles=0
	num_true_chr_clustered_bubbles=0
	num_haplo_clustered_bubbles=0
	num_true_haplo_clustered_bubbles=0
	num_haploclust_false_pos=0
	num_haploclust_false_neg=0

	# defining the same values per chromosome
	chrom_num_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_chr_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_garbage_chr_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_true_chr_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_haplo_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_true_haplo_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_false_haplo_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_haploclust_false_pos={chrom:0 for chrom in valid_chroms}
	chrom_num_haploclust_false_neg={chrom:0 for chrom in valid_chroms}
	chrom_switch_haplo={chrom:False for chrom in valid_chroms}
	
	#for bubble_id, bubble in bubbles.items():
	for bubble_id in bubbles:
		bubble = bubbles[bubble_id]

		if bubble.actual_chrom == None or bubble.actual_chrom not in valid_chroms:
			# the bubble does not have a valid actual chrom
			continue

		chrom = bubble.actual_chrom

		chrom_num_bubbles[chrom]+=1

		if bubble.clust == None:
			bubble.pred_type="not_chrom_clust"
			continue

		if bubble.clust not in clust_to_chrom:
			num_garbage_chr_clustered_bubbles+=1
			chrom_num_garbage_chr_clustered_bubbles[chrom]+=1
			bubble.pred_type="garbage_clust"
			continue

		num_chr_clustered_bubbles+=1
		chrom_num_chr_clustered_bubbles[chrom]+=1

		if clust_to_chrom[bubble.clust] == bubble.actual_chrom:
			num_true_chr_clustered_bubbles+=1
			chrom_num_true_chr_clustered_bubbles[chrom]+=1

		if bubble.allele0.pred_haplo == None:
			bubble.pred_type="not_haplo_clust"

			if bubble.allele0.actual_haplo != None:
				bubble.pred_type="haplo_clust_false_neg"
				num_haploclust_false_neg+=1
				chrom_num_haploclust_false_neg[chrom]+=1

			continue

		num_haplo_clustered_bubbles+=1
		chrom_num_haplo_clustered_bubbles[chrom]+=1

		if bubble.allele0.actual_haplo==None:
			bubble.pred_type="haploclust_false_pos"
			num_haploclust_false_pos+=1
			chrom_num_haploclust_false_pos[chrom]+=1
			continue

		elif bubble.allele0.pred_haplo == bubble.allele0.actual_haplo:
			bubble.pred_type="true_haplo_clust"
			chrom_num_true_haplo_clustered_bubbles[chrom]+=1
			continue

		else:
			bubble.pred_type="false_haplo_clust"
			chrom_num_false_haplo_clustered_bubbles[chrom]+=1

	# revise chrom num true haplotypes after switching the haplotypes in chroms if necessary
	for chrom in valid_chroms:
		true_haplo_clust=chrom_num_true_haplo_clustered_bubbles[chrom]
		false_haplo_clust=chrom_num_false_haplo_clustered_bubbles[chrom]
		
		if true_haplo_clust < false_haplo_clust:
			# switch haplotypes in the chrom
			chrom_switch_haplo[chrom]=True
			true_haplo_clust, false_haplo_clust = false_haplo_clust, true_haplo_clust
			num_true_haplo_clustered_bubbles+=true_haplo_clust
			chrom_num_true_haplo_clustered_bubbles[chrom]=true_haplo_clust
			chrom_num_false_haplo_clustered_bubbles[chrom]=false_haplo_clust

	# revise the type of the bubble if the haplotype is switched in the chromosome
	for bubble_id, bubble in bubbles.items():
		if bubble.pred_type=="true_haplo_clust" or bubble.pred_type=="false_haplo_clust":
			if chrom_switch_haplo[bubble.actual_chrom]:
				bubble.pred_type="true_haplo_clust" if bubble.pred_type=="false_haplo_clust" else "false_haplo_clust"

	# writing the performance statistics in the outout file
	with open(output_file, 'w') as out:
		print('*** Note: false positive haplotype clustered bubbles are also counted in haplotype clustering accuracy', file=out)
		print('total number of bubbles =', num_bubbles, file=out)
		print('number of chrom clustered bubbles =', num_chr_clustered_bubbles, ', (', num_chr_clustered_bubbles*100/num_bubbles, ' % of #bubbles)', file=out)
		print('chrom clustering accuracy =', num_true_chr_clustered_bubbles*100/num_chr_clustered_bubbles, file=out)
		print('number of haplo clustered bubbles = ', num_haplo_clustered_bubbles, ', (', num_haplo_clustered_bubbles*100/num_bubbles, '% of #bubbles)', file=out)
		print('haplo clustering accuracy =', num_true_haplo_clustered_bubbles*100/num_haplo_clustered_bubbles, file=out)
		print('haplo clustering false positive rate =', num_haploclust_false_pos*100/num_bubbles, file=out)
		print('haplo clustering false negative rate =', num_haploclust_false_neg*100/num_bubbles, file=out)


		print('chromosome wise clustering accuracy:')

		print('chrom\t#bubbles\t\
		#chr_clustered_bubbles\tfraction_chr_clustered_bubbles\tchr_clustering_accuracy\t\
		#haplo_clustered_bubbles\tfraction_haplo_clustered_bubbles\thaplo_clustering_accuracy\t\
		false_pos\tfalse_neg', file=out)

		for chrom in valid_chroms:
			print(chrom, '\t', chrom_num_bubbles[chrom], '\t', 
			chrom_num_chr_clustered_bubbles[chrom], '\t', 
			chrom_num_chr_clustered_bubbles[chrom]*100/chrom_num_bubbles[chrom], '\t', 
			chrom_num_true_chr_clustered_bubbles[chrom]*100/chrom_num_chr_clustered_bubbles[chrom], '\t',
			chrom_num_haplo_clustered_bubbles[chrom], '\t', 
			chrom_num_haplo_clustered_bubbles[chrom]*100/chrom_num_bubbles[chrom], '\t', 
			chrom_num_true_haplo_clustered_bubbles[chrom]*100/chrom_num_haplo_clustered_bubbles[chrom], '\t', 
			chrom_num_haploclust_false_pos[chrom]*100/chrom_num_bubbles[chrom], '\t', 
			chrom_num_haploclust_false_neg[chrom]*100/chrom_num_bubbles[chrom], '\t', 		
			file=out)

	print('elapsed time =', time.time()-start_time)	


def output_bubbles_haplo_dist(bubbles, output_file, with_km=True):

# TODO: assert that alignments are present
	'''
	'''

	start_time = time.time()
	print('outputting the bubbles haplo edit dist')

	with open(output_file, 'w') as out:
		km = "km" if with_km else ""
		print('chrom\tbubble_id\tdist_h1\tdist_h2\tnum_aln_reads_h1\tnum_aln_reads_h2\tactual_type\tpred_type\t'+str(km), file=out)

		for bubble_id, bubble in bubbles.items():

			h0_edit_dist, h1_edit_dist = bubble.allele0.get_haplotypes_edit_dist() #bubble.get_haplotypes_edit_dist()			
			h0_num_aln_reads, h1_num_aln_reads = bubble.get_haplo_num_aligned_reads()

			km = bubble.allele0.km+bubble.allele1.km if with_km else ""

			print(str(bubble.actual_chrom), '\t', bubble_id, '\t', 
						h0_edit_dist, '\t', h1_edit_dist, '\t', 
						h0_num_aln_reads, '\t', h1_num_aln_reads, '\t', 
						bubble.actual_type, '\t', bubble.pred_type, '\t', km, file=out)

	print('elapsed time =', time.time()-start_time)


def evaluate_long_read_clustering(long_reads, output_file):

	pass
	'''
	'''

	start_time = time.time()
	print('evaluating long reads clustering')

	num_long_reads=len(long_reads)
	num_haplo_clustered_long_reads=0
	num_true_haplo_clustered_long_reads=0
	num_haploclust_false_pos=0
	num_haploclust_false_neg=0

	# defining the same values per chromosome
	chrom_num_long_reads={chrom:0 for chrom in valid_chroms}
	chrom_num_haplo_clustered_long_reads={chrom:0 for chrom in valid_chroms}
	chrom_num_true_haplo_clustered_long_reads={chrom:0 for chrom in valid_chroms}
	chrom_num_false_haplo_clustered_long_reads={chrom:0 for chrom in valid_chroms}
	chrom_num_haploclust_false_pos={chrom:0 for chrom in valid_chroms}
	chrom_num_haploclust_false_neg={chrom:0 for chrom in valid_chroms}
	chrom_switch_haplo={chrom:False for chrom in valid_chroms}
	
	for read_name, long_read in long_reads.items():

		if long_read.actual_chrom == None or long_read.actual_chrom not in valid_chroms:
			# the long_read does not have a valid actual chrom
			continue

		chrom = long_read.actual_chrom

		chrom_num_long_reads[chrom]+=1

		if long_read.pred_haplo == None:
			long_read.pred_type="not_haplo_clust"

			if long_read.actual_haplo != None:
				long_read.pred_type="haplo_clust_false_neg"
				num_haploclust_false_neg+=1
				chrom_num_haploclust_false_neg[chrom]+=1

			continue

		num_haplo_clustered_long_reads+=1
		chrom_num_haplo_clustered_long_reads[chrom]+=1

		if long_read.actual_haplo==None:
			long_read.pred_type="haploclust_false_pos"
			num_haploclust_false_pos+=1
			chrom_num_haploclust_false_pos[chrom]+=1
			continue

		elif long_read.pred_haplo == long_read.actual_haplo:
			long_read.pred_type="true_haplo_clust"
			chrom_num_true_haplo_clustered_long_reads[chrom]+=1
			continue

		else:
			long_read.pred_type="false_haplo_clust"
			chrom_num_false_haplo_clustered_long_reads[chrom]+=1

	# revise chrom num true haplotypes after switching the haplotypes in chroms if necessary
	for chrom in valid_chroms:
		true_haplo_clust=chrom_num_true_haplo_clustered_long_reads[chrom]
		false_haplo_clust=chrom_num_false_haplo_clustered_long_reads[chrom]
		
		if true_haplo_clust < false_haplo_clust:
			# switch haplotypes in the chrom
			chrom_switch_haplo[chrom]=True
			true_haplo_clust, false_haplo_clust = false_haplo_clust, true_haplo_clust
			num_true_haplo_clustered_long_reads+=true_haplo_clust
			chrom_num_true_haplo_clustered_long_reads[chrom]=true_haplo_clust
			chrom_num_false_haplo_clustered_long_reads[chrom]=false_haplo_clust

	# revise the type of the long_read if the haplotype is switched in the chromosome
	for read_name, long_read in long_reads.items():
		if long_read.pred_type=="true_haplo_clust" or long_read.pred_type=="false_haplo_clust":
			if chrom_switch_haplo[long_read.actual_chrom]:
				long_read.pred_type="true_haplo_clust" if long_read.pred_type=="false_haplo_clust" else "false_haplo_clust"

	# writing the performance statistics in the outout file
	with open(output_file, 'w') as out:
		print('*** Note: false positive haplotype clustered long_reads are also counted in haplotype clustering accuracy', file=out)
		print('total number of long_reads =', num_long_reads, file=out)
		#print('number of chrom clustered long_reads =', num_chr_clustered_long_reads, ', (', num_chr_clustered_bubbles*100/num_long_reads, ' % of #long_reads)', file=out)
		#print('chrom clustering accuracy =', num_true_chr_clustered_long_reads*100/num_chr_clustered_long_reads, file=out)
		print('number of haplo clustered long_reads = ', num_haplo_clustered_long_reads, ', (', num_haplo_clustered_long_reads*100/num_long_reads, '% of #long_reads)', file=out)
		
		print('haplo clustering accuracy =', num_true_haplo_clustered_long_reads*100/num_haplo_clustered_long_reads, file=out)
		print('haplo clustering false positive rate =', num_haploclust_false_pos*100/num_long_reads, file=out)
		print('haplo clustering false negative rate =', num_haploclust_false_neg*100/num_long_reads, file=out)


		print('chromosome wise clustering accuracy:')

		print('chrom\t#long_reads\t\
		#haplo_clustered_long_reads\tfraction_haplo_clustered_long_reads\thaplo_clustering_accuracy\t\
		false_pos\tfalse_neg', file=out)

		for chrom in valid_chroms:
			print(chrom, '\t', chrom_num_long_reads[chrom], '\t', 
			#chrom_num_chr_clustered_long_reads[chrom], '\t', 
			#chrom_num_chr_clustered_long_reads[chrom]*100/chrom_num_long_reads[chrom], '\t', 
			#chrom_num_true_chr_clustered_long_reads[chrom]*100/chrom_num_chr_clustered_long_reads[chrom], '\t',
			chrom_num_haplo_clustered_long_reads[chrom], '\t', 
			chrom_num_haplo_clustered_long_reads[chrom]*100/chrom_num_long_reads[chrom], '\t', 
			chrom_num_true_haplo_clustered_long_reads[chrom]*100/chrom_num_haplo_clustered_long_reads[chrom], '\t', 
			chrom_num_haploclust_false_pos[chrom]*100/chrom_num_long_reads[chrom], '\t', 
			chrom_num_haploclust_false_neg[chrom]*100/chrom_num_long_reads[chrom], '\t', 		
			file=out)

	print('elapsed time =', time.time()-start_time)	


def output_long_reads_haplo_dist(long_reads, output_file):

# TODO: assert that alignments are present
	'''
	'''

	start_time = time.time()
	print('outputting the long reads haplo edit dist')

	with open(output_file, 'w') as out:
		print('chrom\tread_name\tdist_h1\tdist_h2\tnum_aln_bubbles\tactual_type\tpred_type', file=out)

		for read_name, long_read in long_reads.items():

			# TODO: check whether it is equal to the previous haplo_dist values ...
			long_read.set_haplotypes_edit_dist() #bubble.get_haplotypes_edit_dist()			
			num_aln_reads = len(long_read.alignments)

			print(str(long_read.actual_chrom), '\t', read_name, '\t', 
						long_read.haplo0_edit_dist, '\t', long_read.haplo1_edit_dist, '\t', 
						num_aln_reads, '\t', 
						long_read.actual_type, '\t', long_read.pred_type, file=out)

	print('elapsed time =', time.time()-start_time)


if __name__ == "__main__":
	parser = ArgumentParser(description=__doc__)
	parser.add_argument("--bubble_bam_file", type=str, help="Bubble haplotagged bam file", required=True)
	parser.add_argument("--bubble_phase_file", type=str, help="Bubbles haplo clustering file", required=True)
	parser.add_argument("--bubble_clust_file", type=str, help="Bubbles cluster file", required=True)
	parser.add_argument("--clust_to_chrom_file", type=str, help="Cluster to chrom mapping file", required=True)
	parser.add_argument("--het_kmers_files", nargs='*', help="The set of bubble/long reads het kmers files", required=True)
	parser.add_argument("--long_read_haplotagged_bam_files", nargs='*', help="The set of long reads haplotagged bam files", required=True)
	parser.add_argument("--long_read_phase_files", nargs='*', help="The set of long reads haplo clustering files", required=True)
	parser.add_argument("--bubbles_clusering_evaluation_file", type=str, help="The output bubbles clustring evaluation file", required=True)
	parser.add_argument("--bubbles_haplo_edit_dist_file", type=str, help="The output bubbles' haplo edit dist file", required=True)
	parser.add_argument("--long_read_phase_evaluation_file", type=str, help="The output long reads phasing evaluation file", required=True)
	parser.add_argument("--long_reads_haplo_edit_dist_file", type=str, help="The output bubbles' haplo edit dist file", required=True)
	parser.add_argument("--with_km", action='store_true', help="True if bubbles have km information in their name")
	args = parser.parse_args()

	#counts(args.sample, args.input_bam, args.input_bed, args.counts_output)
	with_km = True if "with_km" in args else False
	print('with_km', 	with_km)

	bubbles = get_bubbles(args.bubble_bam_file, with_km)
	clust_to_chrom = get_clust_to_chrom(args.clust_to_chrom_file)
	add_bubble_clust(args.bubble_clust_file, bubbles)
	add_bubble_allele_pred_haplo(args.bubble_phase_file, bubbles)
	long_reads = get_long_reads(args.long_read_haplotagged_bam_files)
	add_long_reads_pred_haplotype(args.long_read_phase_files, long_reads)
	set_alignments(args.het_kmers_files, bubbles, long_reads)
	evaluate_bubble_clustering(bubbles, clust_to_chrom, args.bubbles_clusering_evaluation_file)
	output_bubbles_haplo_dist(bubbles, args.bubbles_haplo_edit_dist_file, with_km)
	evaluate_long_read_clustering(long_reads, args.long_read_phase_evaluation_file)
	output_long_reads_haplo_dist(long_reads, args.long_reads_haplo_edit_dist_file)

