from bubble_long_read_alignment import *
import pysam
import time
import pdb

global valid_chroms
valid_chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX']


def get_bubbles(bubble_haplotagged_bam_file, with_km=True):

	'''
	Reads bubbles happlotagged bam files and returns the list of bubbles

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
				assert (read_name in phased_long_reads), 'long read ' + read_name + ' is not present in long reads'

				bubble = bubbles[bubble_id]
				long_read = phased_long_reads[read_name]

				bubble_allele = bubble.allele0 if bubble_al==0 else bubble.allele1

				assert (bubble_allele != None), 'bubble ' + str(bubble_id) + ' allele ' + str(bubble_al) + ' is None'
				
				aln = Alignment(long_read=long_read, bubble_allele=bubble_allele, bubble_kmer=bubble_kmer, long_read_kmer=long_read_kmer, edit_dist=edit_dist)

	print('elapsed time =', time.time()-start_time)	