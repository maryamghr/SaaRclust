from __future__ import division
from bubble_long_read_alignment import *
import pysam
import numpy
import time
import pdb
from argparse import ArgumentParser

global valid_chroms
valid_chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX']


## supplementary functions
def print_dict_head(dic, num):
	for i in range(num):
		key = list(dic.keys())[i]
		print(key, ':', dic[key])

def print_dp_table(dp_table, str1=None, str2=None, max_dist_digits=2):
	def adjust_space(char, max_digits):
		n_digits = len(str(char))
		add_space = " "*(max_digits - n_digits)
		return add_space+str(char)
	
	if str2 != None:
		print("-", [adjust_space(char, max_dist_digits) for char in '-'+str2])
		
	for i in range(len(dp_table)):
		char = ""
		if str1 != None:
			char = "-" if i==0 else str1[i-1]			
			
		print(char, [adjust_space(dist, max_dist_digits) for dist in dp_table[i]])
		
def needleman_wunsch(str1, str2):
	m, n = len(str1), len(str2)
	edit_dist = [[j for j in range(n+1)]]
	
	for i in range(1, m+1):
		edit_dist.append([i])
		
		for j in range(1, n+1):				
			d = 0 if str1[i-1]==str2[j-1] else 1
			
			edit_dist[i].append(min(d+edit_dist[i-1][j-1], \
											1+edit_dist[i-1][j], # str1[i] aligned with gap \
											1+edit_dist[i][j-1]
											))
	assert (len(edit_dist)==m+1), 'number of rows in edit distance should be ' + str(m+1)
	for i in range(m+1):
		assert (len(edit_dist[i])==n+1), 'number of columns in edit distance should be ' + str(n+1)
		
	# backtracking
	str1_aln, aln, str2_aln = "", "", ""
	
	i, j = m, n
	
	while i > 0 or j > 0:
		#print('i, j =', i, j)
		#print(str1_aln+"\n"+aln+"\n"+str2_aln)
	
		back_dist = []
		if j != 0:
			back_dist.append(1+edit_dist[i][j-1])
			
		if i != 0:
			back_dist.append(1+edit_dist[i-1][j])
		
		if i != 0 and j != 0:
			d = 0 if str1[i-1]==str2[j-1] else 1
			back_dist.append(d+edit_dist[i-1][j-1])
			
		min_back_dist = min(back_dist)
		
		if i != 0 and j != 0 and min_back_dist == d+edit_dist[i-1][j-1]:
			# pos i in str1 is aligned with pos j in str2
			str1_aln = str1[i-1] + str1_aln
			str2_aln = str2[j-1] + str2_aln
			aln_char = "|" if str1[i-1]==str2[j-1] else "."
			aln = aln_char + aln
			i, j = i-1, j-1
			continue
		
		elif j != 0 and min_back_dist == 1+edit_dist[i][j-1]:
			# pos j in str2 is aligned with gap
			str1_aln = "-" + str1_aln
			str2_aln = str2[j-1] + str2_aln
			aln = " " + aln
			j = j-1
			continue
			
		else: # i!=0 and min_back_dist == 1+edit_dist[i-1][j]:
			# pos i in str1 is aligned with gap
			str1_aln = str1[i-1] + str1_aln
			str2_aln = "-" + str2_aln
			aln = " " + aln
			i = i-1
			
	print(str1_aln)
	print(aln)
	print(str2_aln)

	return edit_dist[m][n]
	
	
def get_alignment_from_cigar(ref, query, aln_ref_start_pos, aln_query_start_pos, cigar):
	'''
	This function returns the alignment from a cigar string
	
	Args:
		ref: the reference sequence
		query: the query sequence
		aln_ref_start_pos: 0-based alignment start position in the reference sequence
		aln_query_start_pos: 0-based alignment start position in the query sequence
		cigar: A string describing how the query sequence aligns to the reference sequence. It has integers followed by characters \in "MIDNSHP=X".
			The characters have the following meaning:
				M: alignment match (consumes_query=yes, consumes_reference=yes)
				I: insertion to the reference (consumes_query=yes, consumes_reference=no)
				D: deletion from the reference (consumes_query=no, consumes_reference=yes)
				N: skipped region from the reference (consumes_query=no, consumes_reference=yes)
				S: soft clipping (consumes_query=yes, consumes_reference=no)
				H: hard clipping (consumes_query=no, consumes_reference=no)
				P: padding: silent deletion from padded reference (consumes_query=no, consumes_reference=no)
				=: sequence match (consumes_query=yes, consumes_reference=yes)
				X: sequence mismatch (consumes_query=yes, consumes_reference=yes)
					
			consumes_query and consumes_reference indicate whether the CIGAR operation causes the alignment to step along the query sequence and the reference sequence respectively
				
	Returns:
		a tuple containing aligned reference, alignment sequence, and aligned query
	'''

	cigar_operations = 'MIDNSHP=X'	
	
	# defining consumes_query and consumes_ref for all cigar operations
	
	consumes_query = {'M': True, 'I': True , 'D': False, 'N': False, 'S': True , 'H': False, 'P': False, '=': True, 'X': True}
	consumes_ref =   {'M': True, 'I': False, 'D': True , 'N': True , 'S': False, 'H': False, 'P': False, '=': True, 'X': True}
		
	ref_pos = aln_ref_start_pos
	query_pos = aln_query_start_pos
	
	aln = ''
	ref_aln = ''
	query_aln = ''
		
	length = ''
	for i in range(len(cigar)):
		if cigar[i].isdigit():
			length = length + cigar[i]
			
		else:
			assert(len(length) > 0), 'there should not be two cigar operation characters next to each other'
			assert(cigar[i] in cigar_operations), "cigar " + cigar + " does not have valid characters"
			
			length = int(length)
			
			aln_char = '|'
			
			if consumes_ref[cigar[i]]:
				ref_aln += ref[ref_pos : ref_pos+length]
				ref_pos += length
			else:
				ref_aln += "-"*length
				aln_char = " "
				
			if consumes_query[cigar[i]]:
				query_aln += query[query_pos : query_pos+length]
				query_pos += length
			else:
				query_aln += "-"*length
				aln_char = " "
				
			aln = aln + aln_char * length

	assert(len(ref_aln)==len(query_aln)), 'the lengths of the aligned sequences should be the same'
	
	return ref_aln, aln, query_aln
	
	
def get_reference_aln_substr(ref, query, aln_ref_start_pos, aln_query_start_pos, cigar, query_start, query_end):
	'''
	This function returns the reference subsequence aligned to the given interval of the query sequence (query start and end are both included in the interval)
	
	Args:
		ref: 						@inherited from get_alignment_from_cigar
		query: 					@inherited from get_alignment_from_cigar
		aln_ref_start_pos: 	@inherited from get_alignment_from_cigar
		aln_query_start_pos: @inherited from get_alignment_from_cigar
		cigar: 					@inherited from get_alignment_from_cigar
		query_start: 0-based start of the interval of interest in the query sequence
		query_end  : 0-based end   of the interval of interest in the query sequence
	
	Returns:
		the reference subsequence aligned to the given interval of the query sequence
	'''

	# getting the alignments from the cigar string
	
	ref_aln, aln, query_aln = get_alignment_from_cigar(ref, query, aln_ref_start_pos, aln_query_start_pos, cigar)
	
	# computing the reference substr of interest
	
	ref_substr = ""
	query_pos = 0
	for i in range(len(ref_aln)):
		if query_start <= query_pos <= query_end and ref_aln[i] != "-":
			ref_substr += ref_aln[i]
		
		if query_aln[i] != "-":
			# we read one char from query
			query_pos += 1
			
	return ref_substr
	

####################################################################################


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

				bubble_id, bubble_al, read_name, bubble_kmer, long_read_kmer, edit_dist = \
				int(sp[0]), int(sp[1]), sp[2], sp[3], sp[4], int(sp[5])

				read_name_sp = read_name.split('/ccs')
				read_name = read_name_sp[0]+'/ccs'

				assert (bubble_id in bubbles), 'bubble ' + str(bubble_id) + ' is not present in the bubbles'
				assert (read_name in long_reads), 'long read ' + read_name + ' is not present in long reads'

				bubble = bubbles[bubble_id]
				long_read = phased_long_reads[read_name]

				bubble_allele = bubble.allele0 if bubble_al==0 else bubble.allele1

				assert (bubble_allele != None), 'bubble ' + str(bubble_id) + ' allele ' + str(bubble_al) + ' is None'
				
				aln = Alignment(long_read, bubble_allele, bubble_kmer, long_read_kmer, edit_dist)

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
	

def output_sampled_long_reads(num_sample, edit_dist_fraction_range, file_name):
	n = 0
	
	with open(file_name, 'w') as out:
		print("bubbleName\tbubbleAllele\tPBname\tbubbleKmer\tPBkmer\tkmersEditDistance\tbubble_alle_pred_haplo\tlong_read_pred_haplo", file=out)
		
		for read_name, long_read in long_reads.items():
			if n > num_sample:
					break
					
			d0, d1 = long_read.haplo0_edit_dist, long_read.haplo1_edit_dist
			if d0+d1 == 0:
				continue
			
			dist_frac = min(d0, d1)/(d0+d1)
		
			if edit_dist_fraction_range[0] < dist_frac < edit_dist_fraction_range[1]:
				n += 1
				
				for aln in long_read.alignments:
						print(aln.output_print(), file=out)
						
				print('h_dist0 =', d0, file=out)
				print('h_dist1 =', d1, '\n', file=out)
						
			

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
	parser.add_argument("--bubbles_haplo_edit_dist_file", type=str, help="The output bubbles haplo edit dist file", required=True)
	parser.add_argument("--long_read_phase_evaluation_file", type=str, help="The output long reads phasing evaluation file", required=True)
	parser.add_argument("--long_reads_haplo_edit_dist_file", type=str, help="The output long reads haplo edit dist file", required=True)

	# testing output files for observation 
	parser.add_argument("--long_reads_with_small_frac_haplo_edit_dist", type=str, help="The output sample long reads with small fraction of haplo_edit_dist", required=True)
	parser.add_argument("--long_reads_with_peak_frac_haplo_edit_dist", type=str, help="The output sample long reads with peak fraction of haplo_edit_dist", required=True)
	######################################
	
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
	
	# sampling long reads:
	output_sampled_long_reads(100, (0.35, 0.4), args.long_reads_with_small_frac_haplo_edit_dist)
	output_sampled_long_reads(10, (0, 0.1), args.long_reads_with_peak_frac_haplo_edit_dist)

