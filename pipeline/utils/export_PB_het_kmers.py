import sys
import copy
import gzip



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
	
	assert(len(het_pos) > 0), "bubble " + str(bubble_id) + ":there should be at least one heterozygoud site."
	assert(het_pos[0] >= q), "bubble " + str(bubble_id) + ":the heterozygous position should be larger than q."
	# FIXME: the following assert fails for bubble 1154957. Fix this problem and remove the following if:
	if het_pos[-1] >= len(chain0_seq)-q:
		print("Warning in bubble " + str(bubble_id) + ":the heterozygous position should be larger than length of the bubble chains minus q.")
		return

	assert(het_pos[-1] < len(chain0_seq)-q), "bubble " + str(bubble_id) + ":the heterozygous position should be larger than length of the bubble chains minus q."

	bubble_het_positions[bubble_id] = het_pos
	bubble_allele_to_kmers[(bubble_id, 0)] = [chain0_seq[i-q:(i+q+1)] for i in het_pos]
	bubble_allele_to_kmers[(bubble_id, 1)] = [chain1_seq[i-q:(i+q+1)] for i in het_pos]
	

def find_query_interval(ref_start, ref_end, alignment_start_pos, cigar):
	'''
	This function returns the index of the query sequence aligned to the ref_index in the reference sequence
	
	Args:
		ref_start: 0-based start index in the referece sequence
		ref_end  : 0-based end   index in the referece sequence
		alignment_start_pos: 0-based alignment start position in the reference sequence
		cigar: A string describing how the query sequence maps to the reference sequence. It has integers followed by characters \in "MIDNSHP=X".
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
	
	# remove the beginning part
	#cigar = cigar.split('cg:Z:')[1]
	
	ref_pos = alignment_start_pos
	query_pos = 0
	
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
	
	query_start_index, query_end_index = 0, 0
	
	for i in range(len(aligned_ref_positions)):
		if aligned_ref_positions[i] == ref_start:
			query_start_index = i
		if aligned_ref_positions[i] == ref_end:
			query_end_index = i
			break
			
	# move the query start index to the right as long as we have gap in the query
	while (query_start_index < len(aligned_query_positions) and aligned_query_positions[query_start_index]=='-'):
		query_start_index += 1
	
	# move the query end index to the left as long as we have gap in the query
	while (query_start_index >= 0 and aligned_query_positions[query_end_index]=='-'):
		query_end_index -= 1		
	
	#return [aligned_ref_positions, aligned_query_positions]	
	return (aligned_query_positions[query_start_index], aligned_query_positions[query_end_index])
	
	
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
	
	
def read_strand_states(strand_state_files):
	'''
	Reads phased strand states from the input file
	
	Args:
		strand_state_file: the file containing phased strand states 
	'''
	
	lib_clust_to_haplo = {}

	for f in strand_state_files:
		file_name = f.split("/")[-1]
		lib_name=file_name.split("_haplo_strand_states")[0]
		with open(f) as states:
			for line in states:
				# process the line if it is not empty nor the header line
				if line!="" and line[0]=="V":
					sp = line.split()
					lib_clust_to_haplo[(lib_name, sp[0])]=int(sp[1])-1
					
	return lib_clust_to_haplo

	
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
		# FIXME: bubbble 3390 in ccs reads has multiple alignments

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


def update_ss_kmer_interval(ss_kmer_interval, ss_end, kmer, alt_kmer):
	'''
	Updates ss_kmer_interval and kmer and alt_kmer if the ss read does not cover the kmer properly
	'''
	
	if ss_kmer_interval[0] < 0:
		# het pos is at the beginning of the ss read
		# cut the beginning of kmer and alt kmer
		kmer = kmer[-ss_kmer_interval[0]:]
		alt_kmer = alt_kmer[-ss_kmer_interval[0]:]
		ss_kmer_interval[0] = 0
						
	if ss_kmer_interval[1] > ss_end:
		# het pos is at the end of the ss read
		# cut the beginning of kmer and alt kmer
		cut = ss_kmer_interval[1] - ss_end
		kmer = kmer[:-cut]
		alt_kmer = alt_kmer[:-cut]
		ss_kmer_interval[1] = ss_end



def read_SS_bubble_map(ss_bubble_map_files, bubble_het_positions, bubble_allele_to_kmers, q):
	'''
	Reads the ss to bubble mapping information and returns the SS reads heterozygous kmers information
	
	Args:
		ss_bubble_map_files: a list of name of files containing the exact mapping information of SS reads to SNV bubbles
		bubble_het_positions: A dictionary that maps a bubble id to its set of heterozygous positions
		bubble_allele_to_kmers: A dictionary that maps a pair (bubble_id, allele) to its set of heterozygous kmers
		q: the length of the extention to the right and left to read from the heterozygous position (=k/2-1 for a kmer)
	'''
	
	ss_to_kmer_altkmer={}
	ss_to_kmer_pos_and_interval={}
	ss_to_clust={}
	ss_bubble_map_dir={}


	for input_file in ss_bubble_map_files: #snakemake.input["SS_bubble_map"]:
		with open(input_file) as f:
			for line in f:
				if line=="":
					break
				
				# skip the header line
				if line.startswith("SSname"):
					continue
					
				sp = line.split()
				
				ss_name, ss_lib, ss_clust, bubble_id, allele, is_reverse = sp[0], sp[1], sp[2], int(sp[3]), int(sp[4]), sp[7]
				map_dir = -1 if is_reverse=="TRUE" else 1
				
				if (bubble_id, allele) not in bubble_allele_to_kmers:
					# the bubble is not clustered
					continue

				# checking which het pos in the bubble is covered by the ss read
				bubble_start, ss_start, aln_len = int(sp[8])-1, int(sp[9])-1, int(sp[10])
				ss_end = ss_start + aln_len -1 # 0-based
				 

				for h in range(len(bubble_het_positions[bubble_id])):
					het_pos = bubble_het_positions[bubble_id][h]
					if het_pos < bubble_start and het_pos - bubble_start >= aln_len:
						# heterozygous position het_pos is not covered by the strand-seq read
						continue
					
					kmer = bubble_allele_to_kmers[(bubble_id, allele)][h]
					alt_kmer = bubble_allele_to_kmers[(bubble_id, 1-allele)][h]
					
					ss_het_pos = het_pos - bubble_start + ss_start

					if map_dir == -1:
						kmer = reversecomp(kmer)
						alt_kmer = reversecomp(alt_kmer)					
						ss_start, ss_end = ss_len-ss_end-1, ss_len-ss_start-1
						ss_het_pos = ss_len - ss_het_pos - 1
					
					ss_kmer_interval = [ss_het_pos-q, ss_het_pos+q]
					
					update_ss_kmer_interval(ss_kmer_interval, ss_end, kmer, alt_kmer)
					assert(ss_kmer_interval[0]>=ss_start and ss_kmer_interval[1]<=ss_end), 'kmer interval should be in the part of ss read that is mapped to the bubble'

					ss_to_kmer_altkmer[ss_name] = (kmer, alt_kmer)
					ss_to_kmer_pos_and_interval[ss_name] = (ss_het_pos, ss_kmer_interval)
						
				ss_to_clust[ss_name]=ss_clust
				ss_bubble_map_dir[ss_name]=map_dir
				
	return ss_to_kmer_altkmer, ss_to_kmer_pos_and_interval, ss_to_clust#, ss_bubble_map_dir
	
	

# pb_to_kmers is a dictionary that maps each pb read name to a list of two lists: the first list contains h0 kmers with their intervals and the second one contains h1 kmers with their intervals in the pb read

def get_pb_kmers(minimap_file, ss_to_kmer_altkmer, ss_to_kmer_pos_and_interval, ss_to_clust, lib_clust_to_haplo):
	
	pb_to_kmers = {}
	pb_to_kmer_intervals = {}

	print("computing pb kmer intervals...")

	#with gzip.open(snakemake.input["SS_PB_minimap"], 'rb') as minimap:
	with gzip.open(minimap_file, 'rb') as minimap:
		for line in minimap:
			line = line.decode("utf-8")
			if line.strip() == "":
				break
			sp=line.split()
			# skip the header line
			if sp[0] == "PBreadNames":
				continue

			pb_name, ss_name, lib_name, ss_len, ss_start, ss_end = sp[0], sp[1], sp[2], int(sp[6]), int(sp[7]), int(sp[8])-1	
			strand, pb_start, pb_end, aln_len, cigar, pb_clust = sp[9],int(sp[14]), int(sp[15])-1, int(sp[17]), sp[18], sp[19]
			
			# remove the beginning part of the cigar string
			cigar = cigar.split('cg:Z:')[1]
		
			# skip processing this line if the ss read is not clustered
			if ss_name not in ss_to_clust:
				continue

	
			ss_clust=ss_to_clust[ss_name]
			map_dir = 1
			
			if strand == "-":
				map_dir = -1
				# update the ss alignment interval (convert reverse(ss) interval to ss interval)
				ss_start, ss_end = ss_len-1-ss_end, ss_len-1-ss_start

			if (lib_name, ss_clust) not in lib_clust_to_haplo:
				# the SS lib is not WC in the cluster
				continue

			h = lib_clust_to_haplo[(lib_name, ss_clust)]
			if ss_name in ss_to_kmer_altkmer: # ss read covers a heterozygous position
				ss_kmers = copy.copy(ss_to_kmer_altkmer[ss_name])
				ss_het_pos = copy.copy(ss_to_kmer_pos_and_interval[ss_name][0])
				ss_kmer_interval = copy.copy(ss_to_kmer_pos_and_interval[ss_name][1])
				
				#aln_dir = map_dir*copy.copy(ss_bubble_map_dir[ss_name])

		
				#if aln_dir != 1: # either ss to bubble or ss to pb alignment direction is reverse
				#	ss_het_pos=ss_len-ss_het_pos-1
				#	ss_kmer_interval[0] = ss_len - ss_kmer_interval[0]-1
				#	ss_kmer_interval[1] = ss_len - ss_kmer_interval[1]-1
				#	ss_kmer_interval.reverse()

				if ss_het_pos < ss_start or ss_het_pos > ss_end-1:
					# the aligned part of the ss read does not contain the het position
					continue

				# update ss_kmer_interval if it is not fully aligned to PB read
				ss_kmer_interval[0] = max(ss_kmer_interval[0], ss_start)
				ss_kmer_interval[1] = min(ss_kmer_interval[1], ss_end-1)

				if pb_name not in pb_to_kmers:
					pb_to_kmers[pb_name]=[[],[]] # two lists corresponding to h0 and h1 kmers

				# computing the kmer interval in the pb read
		
				pb_kmer_interval = find_reference_interval(ss_kmer_interval[0], ss_kmer_interval[1], pb_start, ss_start, cigar)
				
				if pb_kmer_interval[0] > pb_kmer_interval[1]:
					continue

				if pb_name not in pb_to_kmer_intervals:
					pb_to_kmer_intervals[pb_name]=[pb_kmer_interval]
				else:
					pb_to_kmer_intervals[pb_name].append(pb_kmer_interval)	
					
				if map_dir==1:
					pb_to_kmers[pb_name][h].append(ss_kmers[0])
					pb_to_kmers[pb_name][1-h].append(ss_kmers[1])
					
				else:
					pb_to_kmers[pb_name][h].append(reversecomp(ss_kmers[0]))
					pb_to_kmers[pb_name][1-h].append(reversecomp(ss_kmers[1]))
					
	return pb_to_kmers, pb_to_kmer_intervals
	
	
# reading pb fasta file and getting het kmer subsequences


def output_pb_het_kmers(fasta_file, output_file, pb_to_kmers, pb_to_kmer_intervals):
	#with open(snakemake.input["PB_fasta"]) as f:
		#with open(snakemake.output[0], 'w') as out:
	print("outputting pb kmers...")
	with open(fasta_file) as f:
		with open(output_file, 'w') as out:			
			print("PB_name\tPB_subseq\th1_kmer\th2_kmer", file=out)
			for line in f:
				if line == "":
					break

				if line[0]==">":
					name = line.strip()[1:]
					name = '_'.join(name.split("_")[:-3])
				else:
					seq = line.strip()
					if name in pb_to_kmers:
						# remove the following if
						#if name != 'm64018_190216_000235/50595562/ccs':
						#	continue
						#if name == 'm64018_190216_000235/50595562/ccs':
						#	print("number of kmers=", len(pb_to_kmers[name][0]))
							#print("pb_kmer_interval=", pb_kmer_interval)
						for j in range(len(pb_to_kmers[name][0])):
							pb_kmer_interval = pb_to_kmer_intervals[name][j]
							subseq=seq[pb_kmer_interval[0]:(pb_kmer_interval[1]+1)]

							#print("pb_to_kmers=", pb_to_kmers[name])
							#print("pb_kmer_interval=", pb_kmer_interval)
						

							#if name == 'm64018_190216_000235/50595562/ccs':
							#	print("subseq=", subseq, "h1_kmer=", pb_to_kmers[name][0][j][0], "h2_kmer=", pb_to_kmers[name][1][j][0])
							#	if subseq == "":
							#		print('j=', j)
							#		sys.exit('PB subseq is empty!!!!')
								

							print(name + "\t" + subseq + "\t" + pb_to_kmers[name][0][j] + "\t" + pb_to_kmers[name][1][j], file=out)
 

