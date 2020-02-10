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
				if line!="" and line[0]!="V":
					sp = line.split()
					lib_clust_to_haplo[(lib_name, sp[0])]=int(sp[1])
					
	return lib_clust_to_haplo


def read_strandphaser_strand_states(strand_states_file):
	'''
	Reads phased strand states from the input file
	
	Args:
		strand_states_file: the file containing phased strand states 
	'''
	
	lib_clust_to_haplo = {}

	with open(strand_states_file) as states:
		#skip the header line
		next(states)
		for line in states:
			if line!="":
				sp = line.split()
				lib_clust_to_haplo[(sp[0], sp[1])]=int(sp[2])
	return lib_clust_to_haplo

def get_ss_clust(ss_clust_file):
	ss_to_clust = {}
	with open(ss_clust_file) as f:
		#skip the header line
		next(f)

		for line in f:
			if line == "":
				continue

			ss_clust, ss_name = line.split()
			ss_to_clust[ss_name] = ss_clust

	return ss_to_clust


def phase_bubbles(ss_bubble_map_files, ss_to_clust, lib_clust_to_haplo, bubble_phase_file):
	'''
	Computes the bubble allele corresponding to the first haplotype (h0) and writes it in the output file
	
	Args:
		ss_bubble_map_files: a list of name of files containing the exact mapping information of SS reads to SNV bubbles
		lib_clust_to_haplo:  A dictionaty that maps a pair (ss_lib, clust) to a haplotype (only for wc ss libs)		
		bubble_phase_file: The output file containing two columns for the bubble name and the bubble allele corresponding to the first haplotype, respectively
	'''
	
	bubble_haplo_allele_count = {}
	bubbles_with_invalid_alleles = {}


	for input_file in ss_bubble_map_files:
		with open(input_file) as f:
			for line in f:
				if line=="":
					break
				
				# skip the header line
				if line.startswith("SSname"):
					continue
					
				sp = line.split()
				
				ss_name, ss_lib, bubble_id, bubble_allele, is_reverse = sp[0], sp[1], sp[2], sp[3], sp[4]
				
				if ss_name not in ss_to_clust:
					continue

				ss_clust = ss_to_clust[ss_name]

				if (ss_lib, ss_clust) not in lib_clust_to_haplo:
					continue

				if bubble_allele!='0' and bubble_allele!='1':
					print('Warning: bubble', bubble_id, 'has non binary allele!')
					bubbles_with_invalid_alleles[bubble_id] = True
					continue

				haplo = lib_clust_to_haplo[(ss_lib, ss_clust)]

				haplo_allele = str(haplo) + bubble_allele
				
				if bubble_id not in bubble_haplo_allele_count:
					bubble_haplo_allele_count[bubble_id] = [0]*4

				bubble_haplo_allele_count[bubble_id][int(haplo_allele, 2)] += 1
				

	with open(bubble_phase_file, 'w') as out:
		print("bubbleName\thaplotype0Allele", file=out)
		for bubble_id in bubble_haplo_allele_count:
			if bubble_id in bubbles_with_invalid_alleles:
				continue

			haplo_allele_count = bubble_haplo_allele_count[bubble_id]
			h0_a0 = haplo_allele_count[0]+haplo_allele_count[3] # count(h=0,a=0) + count(h=1,a=1)
			h0_a1 = haplo_allele_count[1]+haplo_allele_count[2]	# count(h=0,a=1) + count(h=1,a=0)

			if h0_a0 == h0_a1:
				continue

			h0_allele = '0' if h0_a0 > h0_a1 else '1'

			print(bubble_id + "\t" + h0_allele, file=out)

