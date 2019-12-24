import sys
import copy
import gzip
from whatshap.align import edit_distance

def reversecomp(seq):
	'''
	returns the reverse complement of string seq
	'''
	
	revcomp = {'a':'t', 'c':'g', 'g':'c', 't':'a', 'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	rc = ''
	for i in range(len(seq)):
		rc = revcomp[seq[i]] + rc
	return rc


def print_dict_head(dic, num):
	for i in range(num):
		key = list(dic.keys())[i]
		print(key, ':', dic[key])


def add_bubble_kmers(bubble_het_positions, bubble_allele_to_kmers, bubble_id, chain0_seq, chain1_seq, q):
	''' adds the set of heterozygous positions and kmers of a bubble to the corresponding dictionaries
		
		Args:
			bubble_het_positions: A dictionary that maps a bubble id to its set of heterozygous positions
			bubble_allele_to_kmers: A dictionary that maps a pair (bubble_id, allele) to its set of heterozygous kmers
			bubble_id: The id of a bubble
			chain0_seq: The sequence of the first bubble chain
			chain1_seq: The sequence of the second bubble chain
			q: the number of characters that should be read in the right and left direction from a het position (=(k-1)/2)
	'''

	bubble_het_positions[bubble_id]=[]
	for al in range(2):
		bubble_allele_to_kmers[(bubble_id, al)]=[]
	
	assert(len(chain0_seq) == len(chain1_seq)), "Error in bubble " + str(bubble_id) + ": the two chains of the bubble should have the same length for SNV bubbles."		

	het_pos = [i for i in range(len(chain0_seq)) if chain0_seq[i] != chain1_seq[i]]
	
	assert(len(het_pos) > 0), "bubble " + str(bubble_id) + ": there should be at least one heterozygous site."
	assert(het_pos[0] >= q), "bubble " + str(bubble_id) + ": the heterozygous position should be larger than q."
	# FIXME: the following assert fails for bubble 1154957. Fix this problem and remove the following if:
	if het_pos[-1] >= len(chain0_seq)-q:
		print("Warning in bubble " + str(bubble_id) + ": the heterozygous position should be larger than length of the bubble chains minus q.")
		return

	assert(het_pos[-1] < len(chain0_seq)-q), "bubble " + str(bubble_id) + ":the heterozygous position should be larger than length of the bubble chains minus q."

	bubble_het_positions[bubble_id] = het_pos
	bubble_allele_to_kmers[(bubble_id, 0)] = [chain0_seq[i-q:(i+q+1)] for i in het_pos]
	bubble_allele_to_kmers[(bubble_id, 1)] = [chain1_seq[i-q:(i+q+1)] for i in het_pos]



def find_reference_interval(query_start, query_end, aln_ref_start_pos, aln_query_start_pos, cigar):
	'''
	This function returns the index of the query sequence aligned to the ref_index in the reference sequence
	
	Args:
		ref_start: 0-based start index in the referece sequence
		ref_end  : 0-based end   index in the referece sequence
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
		a tuple containing 0-based query start and end positions
	'''

	cigar_operations = 'MIDNSHP=X'	
	
	# defining consumes_query and consumes_ref for all cigar operations
	
	consumes_query = {'M': True, 'I': True , 'D': False, 'N': False, 'S': True , 'H': False, 'P': False, '=': True, 'X': True}
	consumes_ref =   {'M': True, 'I': False, 'D': True , 'N': True , 'S': False, 'H': False, 'P': False, '=': True, 'X': True}
		
	ref_pos = aln_ref_start_pos
	query_pos = aln_query_start_pos
	
	aligned_ref_positions = []
	aligned_query_positions = []
	
	
	length = ''
	for i in range(len(cigar)):
		if cigar[i].isdigit():
			length = length + cigar[i]
			
		else:
			assert(len(length) > 0), 'there should not be two cigar operation characters next to each other'
			assert(cigar[i] in cigar_operations), "cigar " + cigar + " does not have valid characters"
			
			length = int(length)
			
			if not consumes_query[cigar[i]] and not consumes_ref[cigar[i]]:
				continue
			
			elif consumes_query[cigar[i]] and consumes_ref[cigar[i]]:
				# move both query and reference positions forward
				aligned_ref_positions   += range(ref_pos,   ref_pos + length)
				aligned_query_positions += range(query_pos, query_pos + length)
				ref_pos += length
				query_pos += length
				
			elif consumes_query[cigar[i]]:
				# move the query sequence forward and align it with gap in the reference sequence
				aligned_ref_positions   += ['-'] * length
				aligned_query_positions += range(query_pos, query_pos + length)
				query_pos += length
				
			else:
				# move the reference sequence forward and align it with gap in the query sequence
				aligned_ref_positions   += range(ref_pos,   ref_pos + length)
				aligned_query_positions += ['-'] * length
				ref_pos += length
				
			length = ''
				
	assert(len(aligned_ref_positions)==len(aligned_query_positions)), 'the lengths of the aligned sequences should be the same'
	
	ref_start, ref_end = 0, 0
	
	for i in range(len(aligned_query_positions)):
		if aligned_query_positions[i] == query_start:
			ref_start = i
		if aligned_query_positions[i] == query_end:
			ref_end = i
			break
			
	# move the ref start to the right as long as we have gap in the ref
	while (ref_start < len(aligned_ref_positions) and aligned_ref_positions[ref_start]=='-'):
		ref_start += 1
	
	# move the ref end to the left as long as we have gap in the ref
	while (ref_end >= 0 and aligned_ref_positions[ref_end]=='-'):
		ref_end -= 1		
	
	#return [aligned_ref_positions, aligned_query_positions]	
	return (aligned_ref_positions[ref_start], aligned_ref_positions[ref_end])


	
def read_het_snv_bubbles(bubbles_file, q):
	'''
	Reads the bubbles fasta file and returns the set of heterozygous positions and kmers for all valid bubbles that are clustered 
	
	Args:
		bubble_file: a fasta file containing SNV bubble chains
		q: the length of the extention to the right and left to read from the heterozygous position (=k/2-1 for a kmer) 
	'''
	
	def add_bubble_pos_and_kmers(bubble_het_positions, bubble_allele_to_kmers, bubble_id, seq0, seq1, q, clust, num_bubble_chains, bubbles_with_invalid_alleles):
		'''
		If the bubble is valid and the clust is not None,
		adds the set of heterozygous positions and kmers of a bubble to the corresponding dictionaries
		
		Args:
			bubbles_with_invalid_alleles: bubbles with alleles other than 0 and 1
		'''

		if clust!="None" and num_bubble_chains == 2 and bubble_id not in bubbles_with_invalid_alleles:				
			assert(len(seq0)==len(seq1)), "Error in bubble " + str(bubble_id) + ": the two chains of the bubble should have the same length for SNV bubbles."				
			add_bubble_kmers(bubble_het_positions, bubble_allele_to_kmers, prev_bubble_id, seq0, seq1, q)
	
	bubble_het_positions = {}
	bubble_allele_to_kmers = {}

	# FIXME: there are some alleles higher than 2 in the input file. They should not exist!!!
	bubbles_with_invalid_alleles = set()

	with open(bubbles_file) as bubbles:
		prev_bubble_id=-1
		prev_clust = "None"
		num_bubble_chains=0
		seq=["", ""]
		for line in bubbles:
			if line=="":
				break

			if line[0]==">":
				sp_line=line.split()
				
				clust, bubble_name = sp_line[1], sp_line[0][1:]
				sp=bubble_name.split("_")
				bubble_id, allele=int(sp[1]), int(sp[3])-1
				
				if allele > 1:
					bubbles_with_invalid_alleles.add(bubble_id)

				if bubble_id != prev_bubble_id:
					if prev_bubble_id != -1:
						# process the previous bubble
						add_bubble_pos_and_kmers(bubble_het_positions, bubble_allele_to_kmers, prev_bubble_id, seq[0], seq[1], q, prev_clust, num_bubble_chains, bubbles_with_invalid_alleles)
					
					prev_bubble_id = bubble_id
					prev_clust = clust
					num_bubble_chains = 1
					
				else:
					num_bubble_chains += 1							
				
					
			else:
				if allele < 2:
					seq[allele] = line.strip()

		add_bubble_pos_and_kmers(bubble_het_positions, bubble_allele_to_kmers, prev_bubble_id, seq[0], seq[1], q, prev_clust, num_bubble_chains, bubbles_with_invalid_alleles)
		
	return bubble_het_positions, bubble_allele_to_kmers


def get_bubble_haplotype(bubble_phase_file):
	bubble_first_haplo_allele = {}
	with open(bubble_phase_file) as f:
		for line in f:
			if line == "":
				break

			sp = line.split()
			
			if sp[0]=="bubbleName":
				continue

			
			bubble_first_haplo_allele[int(sp[0])] = int(sp[1])

	return bubble_first_haplo_allele


def get_bubble_kmers_with_interval(bubble_het_pos, bubble_aln_start, bubble_aln_end, kmer0, kmer1 , q):
	'''
	Updates bubble_kmer_interval
	'''

	start, end = bubble_het_pos-q, bubble_het_pos+q
	bubble_kmer0, bubble_kmer1 = kmer0, kmer1
	start_cut	= bubble_aln_start - start
	end_cut	  = end - bubble_aln_end

	if start_cut > 0:
		start += start_cut
		bubble_kmer0, bubble_kmer1 = bubble_kmer0[start_cut:], bubble_kmer1[start_cut:]

	if start_cut > 0:
		end -= end_cut
		bubble_kmer0, bubble_kmer1 = bubble_kmer0[:-end_cut], bubble_kmer1[:-end_cut]


	return (start, end, bubble_kmer0, bubble_kmer1)

	
def get_pb_name_to_seq(fasta_file):
	pb_name_to_seq = {}
	with open(fasta_file) as f:
		for line in f:
			if line == "":
				break

			if line[0]==">":
				name = line.strip()[1:]
				#name = '_'.join(name.split("_")[:-3])
					
			else:
				seq = line.strip()
				pb_name_to_seq[name] = seq

	return pb_name_to_seq
				

# pb_to_kmers is a dictionary that maps each pb read name to a list of two lists: the first list contains h0 kmers with their intervals and the second one contains h1 kmers with their intervals in the pb read

def output_bubble_and_pb_kmers(minimap_file, bubble_first_haplo_allele, bubble_het_positions, bubble_allele_to_kmers, pb_name_to_seq, q, output_file):
	
	print("computing pb kmer intervals...")

	with gzip.open(minimap_file, 'rb') as minimap:
		with open(output_file, 'w') as out:
			print("bubbleName\tbubbleAllele\tPBname\tbubbleKmer\tPBkmer\tbubbleAlleleHaplotype\tkmersEditDistance", file=out)

			for line in minimap:
				line = line.decode("utf-8")
				if line.strip() == "":
					break
				sp=line.split()
				
				bubble_name, bubble_len, bubble_start, bubble_end, strand, pb_name, pb_start, pb_end = sp[0], int(sp[1]), int(sp[2]), int(sp[3]), sp[4], sp[5], int(sp[7]), int(sp[8])
				bubble_name_sp = bubble_name.split('_')
				bubble_id, allele = int(bubble_name_sp[1]), int(bubble_name_sp[3])

				if strand=='-':
					bubble_start, bubble_end = bubble_len-1-bubble_end, bubble_len-1-bubble_start

				# process only allele 0 because of symmetry in building the dict
				if allele == 1:
					continue

				if bubble_id not in bubble_het_positions:
					# the bubble is not to be used (because of not being clustered, having invalid allele, ...)
					continue

				# remove the beginning part of the cigar string
				cigar = sp[-1] # assuming that always the last tag is cigar
				cigar = cigar.split('cg:Z:')[1]
		
				het_pos = bubble_het_positions[bubble_id]
				bubble_al0_kmers, bubble_al1_kmers = bubble_allele_to_kmers[(bubble_id, 0)], bubble_allele_to_kmers[(bubble_id, 1)]

				assert(len(het_pos)==len(bubble_al0_kmers) and len(bubble_al0_kmers) == len(bubble_al1_kmers)), "Error in bubble " + str(bubble_id) + ": the two chains of the bubble should have the same number of heterozygous kmers."	

				for i in range(len(het_pos)):
					pos, kmer0, kmer1 = het_pos[i], bubble_al0_kmers[i], bubble_al1_kmers[i]


					if strand=='-':
						pos, kmer0, kmer1 = bubble_len-1-pos, reversecomp(kmer0), reversecomp(kmer1)

					if pos < bubble_start or pos > bubble_end:
						# the heterozygous position is not covered by the long read
						continue

					start, end, bubble_kmer0, bubble_kmer1 = get_bubble_kmers_with_interval(pos, bubble_start, bubble_end, kmer0, kmer1 , q)

					pb_kmer_start, pb_kmer_end = find_reference_interval(start, end, pb_start, bubble_start, cigar)

					if pb_kmer_start > pb_kmer_end:
						# the bubble is aligned with only gap in the long read
						continue

					pb_kmer = pb_name_to_seq[pb_name][pb_kmer_start:(pb_kmer_end+1)]
					
					# note that the haplotype of al0 is 0 iff allele of haplotype 0 is al0 (i.e., bubble_al0_haplo = bubble_first_haplo_allele)
					bubble_al0_haplo = '?' if bubble_id not in bubble_first_haplo_allele else bubble_first_haplo_allele[bubble_id]
					bubble_al1_haplo = '?' if bubble_al0_haplo=='?' else 1-bubble_al0_haplo

					edit_dist_bubble_al0_pb = edit_distance(pb_kmer, bubble_kmer0)
					edit_dist_bubble_al1_pb = edit_distance(pb_kmer, bubble_kmer1)

					print(str(bubble_id) + "\t0\t" + pb_name + "\t" + bubble_kmer0 + "\t" + pb_kmer + "\t" + str(bubble_al0_haplo) + "\t" + str(edit_dist_bubble_al0_pb), file=out)
					print(str(bubble_id) + "\t1\t" + pb_name + "\t" + bubble_kmer1 + "\t" + pb_kmer + "\t" + str(bubble_al1_haplo) + "\t" + str(edit_dist_bubble_al1_pb), file=out)
				

