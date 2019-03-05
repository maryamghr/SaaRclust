import sys
#import re

q = snakemake.params["het_kmer_len"]

print("q = ", q)

revcomp = {'a':'t', 'c':'g', 'g':'c', 't':'a', 'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def reversecomp(seq):
	rc = ''
	for i in range(len(seq)):
		rc = revcomp[seq[i]] + rc
	return rc

def add_bubble_kmers(bubble_allele_to_kmers, bubble_id, seq, q):
	for al in range(2):
		bubble_allele_to_kmers[(bubble_id, al)]={}
	d = []
	for i in range(len(seq[0])):
		if seq[0][i] != seq[1][i]:
			for al in range(2):
				bubble_allele_to_kmers[(bubble_id, al)][i]=seq[al][i-q:(i+q+1)]


def map_index(idx, cigar):
	# counting the number of insertions and deletions that happen before each index in the aligned text	
	num_insertion=[]
	num_deletion= []
	insertion = 0
	deletion = 0
	#print("cigar = ", cigar)
	cigar = cigar.split('cg:Z:')[1]
	l = ''

	for i in range(len(cigar)):
		if cigar[i].isdigit():
			l = l + cigar[i]
		else:
			num = int(l)
			if cigar[i]=='M' or cigar[i]=='X' or cigar[i]=='I':
				num_insertion = num_insertion + [insertion]*num
				num_deletion = num_deletion + [deletion]*num
				if cigar[i]=='I':
					insertion = insertion + num
			elif cigar[i]=='D':
				deletion = deletion + num
			
			else:
				print("cigar=", cigar, " has chars except for {M, X, I, D}")
				return ""

			l = ''

	# map text1 index to text2 index (idx2 = idx1 + |D| - |I|)
	print('num_deletion=', num_deletion)
	print('num_insertion=', num_insertion)
	return idx + num_deletion[idx] - num_insertion[idx]


bubble_allele_to_kmers = {}

# FIXME: there are some alleles higher than 2 in the input file. They should not exist!!!
bubbles_with_invalid_alleles = set()

with open(snakemake.input["bubbles"]) as bubbles:
	prev_bubble_id=-1
	prev_clust = "None"
	num_bubble_rep=0
	seq=["", ""]
	for line in bubbles:
		if line=="":
			break

		if line[0]==">":
			sp_line=line.split()
			clust=sp_line[1]
			bubble_name=sp_line[0][1:]
			sp=bubble_name.split("_")
			bubble_id=int(sp[1])
			allele=int(sp[3])-1
			if allele > 2:
				bubbles_with_invalid_alleles.add(bubble_id)

			if bubble_id != prev_bubble_id:			
				if prev_clust!="None" and num_bubble_rep < 3 and prev_bubble_id not in bubbles_with_invalid_alleles:				
					# process the previous bubble
					if len(seq[0])!=len(seq[1]):
						print("bubble", prev_bubble_id, seq[0], seq[1], ", num_bubble_rep = ", num_bubble_rep)
						sys.exit("Error: the sequence lengths are different!!!")
				
					if prev_bubble_id != -1:
						add_bubble_kmers(bubble_allele_to_kmers, prev_bubble_id, seq, q)		
				
				prev_bubble_id = bubble_id
				prev_clust = clust
				num_bubble_rep = 0
		else:
			if allele < 2:
				seq[allele] = line.strip()
				num_bubble_rep = num_bubble_rep + 1

	add_bubble_kmers(bubble_allele_to_kmers, bubble_id, seq, q)

### test
print("bubble_allele_to_kmers[(64414,0)]=", bubble_allele_to_kmers[(64414,0)])
print("bubble_allele_to_kmers[(64414,1)]=", bubble_allele_to_kmers[(64414,1)])
#print(bubble_allele_to_kmers)


ss_to_kmer_altkmer={}
ss_to_kmer_pos_and_interval={}
ss_to_clust={}
ss_bubble_map_dir={}

#with open(snakemake.input["SS_bubble_map"]) as f:
#	for line in f:
#		if line=="":
#			break
#		
#		if line[0]==">":
#			ss_name=line.split()[1].split("_lib")[0]
#		else:
#			sp=line.split()
#			bubble_id=int(sp[0].split("_")[1])
#			allele=int(sp[0].split("_")[3])-1
#			print("ss_name = ", ss_name, "bubble_id = ", bubble_id, ", allele = ", allele)
#			
#			if (bubble_id, allele) not in bubble_allele_to_kmers:
#				continue
#			# checking which het pos in the bubble in covered by the ss read
#			bubble_start = int(sp[1])-1
#			aln_len = int(sp[3])
#			for het_pos in bubble_allele_to_kmers[(bubble_id, allele)]:
#				if het_pos >= bubble_start and het_pos - bubble_start < aln_len:
#					# heterozygous position het_pos is covered by the strand-seq read
#					ss_to_kmer_altkmer[ss_name]=(bubble_allele_to_kmers[(bubble_id, allele)][het_pos], bubble_allele_to_kmers[(bubble_id, 1-allele)][het_pos])
#					break


# alternative
for input_file in snakemake.input["SS_bubble_map"]:
	with open(input_file) as f:
		for line in f:
			if line=="":
				break
		
			sp = line.split()
			if sp[0]!="SSname":
				ss_name=sp[0]
				ss_lib=sp[1]
				ss_clust=sp[2]
				bubble_id=int(sp[3])
				allele=int(sp[4])
				is_reverse=sp[7]
				map_dir = -1 if is_reverse=="TRUE" else 1
				#print("ss_name = ", ss_name, "bubble_id = ", bubble_id, ", allele = ", allele)
			
				if (bubble_id, allele) not in bubble_allele_to_kmers:
					continue

				# checking which het pos in the bubble in covered by the ss read
				bubble_start = int(sp[8])-1
				ss_start = int(sp[9])-1
				aln_len = int(sp[10])
				for het_pos in bubble_allele_to_kmers[(bubble_id, allele)]:
					if het_pos >= bubble_start and het_pos-bubble_start < aln_len:
						# heterozygous position het_pos is covered by the strand-seq read
						kmer = bubble_allele_to_kmers[(bubble_id, allele)][het_pos]
						alt_kmer = bubble_allele_to_kmers[(bubble_id, 1-allele)][het_pos]
						ss_het_pos = het_pos - bubble_start + ss_start
						ss_kmer_interval = [ss_het_pos-q, ss_het_pos+q]
						# update ss_kmer_interval and kmer and alt_kmer if the ss read does not cover the kmer properly
						if ss_kmer_interval[0] < 0:
							# het pos is at the beginning of the ss read
							# cut the beginning of kmer and alt kmer
							kmer = kmer[-ss_kmer_interval[0]:]
							alt_kmer = alt_kmer[-ss_kmer_interval[0]:]
							ss_kmer_interval[0] = 0
							#ss_kmer_interval = (0, ss_kmer_interval[1])
						if ss_kmer_interval[1] > aln_len + ss_start-1:
							# het pos is at the end of the ss read
							# cut the beginning of kmer and alt kmer
							cut = ss_kmer_interval[1] - aln_len + ss_start + 1
							kmer = kmer[:-cut]
							alt_kmer = alt_kmer[:-cut]
							ss_kmer_interval[1] = aln_len + ss_start-1
							#ss_kmer_interval = [ss_kmer_interval[0], aln_len + ss_start-1]
						ss_to_kmer_altkmer[ss_name]=(kmer, alt_kmer)
						ss_to_kmer_pos_and_interval[ss_name]=[ss_het_pos, ss_kmer_interval]
						ss_to_clust[ss_name]=ss_clust
						ss_bubble_map_dir[ss_name]=map_dir
						break

### test
#print("ss_to_kmer_altkmer:")
#print(ss_to_kmer_altkmer)
print("ss_to_kmer_altkmer['HISEQ:217:HAWJ1ADXX:2:2101:8498:46241']=", ss_to_kmer_altkmer['HISEQ:217:HAWJ1ADXX:2:2101:8498:46241'])

# mapping (lib, clust) pairs to haplotype
lib_clust_to_haplo = {}

for f in snakemake.input["SS_haplo_strand_states"]:
	file_name = f.split("/")[-1]
	lib_name=file_name.split("_haplo_strand_states")[0]
	with open(f) as states:
		for line in states:
			# process the line if it is not empty nor the header line
			if line!="" and line[0]=="V":
				sp = line.split()
				lib_clust_to_haplo[(lib_name, sp[0])]=int(sp[1])-1

#TODO: just for testing, remove later!
lib_clust_to_haplo[('NW150212-IV_92', 'V47')]=0
lib_clust_to_haplo[('NW150212-IV_92', 'V1')] =1


### test
#print("lib_clust_to_haplo:")
#print(lib_clust_to_haplo)


#cluster = snakemake.wildcards["cluster"]
#print("cluster = ", cluster)

pb_to_kmers={}

#with open(snakemake.input["SS_PB_minimap"]) as minimap:
#	for line in minimap:
#		if line.strip() == "":
#			break
#		sp=line.split()
#		pb_name=sp[0]
#		ss_name=sp[1]
#		lib_name=sp[2]
#		if (lib_name, cluster) not in lib_clust_to_haplo:
#			# the SS lib is not WC in the cluster
#			continue
#		h = lib_clust_to_haplo[(lib_name, cluster)]
#		if ss_name in ss_to_kmer_altkmer:
#			ss_kmers = ss_to_kmer_altkmer[ss_name]
#			if pb_name not in pb_to_kmers:
#				pb_to_kmers[pb_name]=[set(),set()] # two lists corresponding to h0 and h1 kmers
#
#			pb_to_kmers[pb_name][h].add(ss_kmers[0])
#			pb_to_kmers[pb_name][1-h].add(ss_kmers[1])
#			
#			# TODO: compute pb interval....


#with open(snakemake.input["SS_PB_minimap"]) as minimap:
#	for line in minimap:
#		if line.strip() == "":
#			break
#		sp=line.split()
#		pb_name=sp[0]
#		ss_name=sp[1]
#		lib_name=sp[2]
#		
#		if (lib_name, cluster) not in lib_clust_to_haplo:
#			# the SS lib is not WC in the cluster
#			continue
#		h = lib_clust_to_haplo[(lib_name, cluster)]
#		if ss_name in ss_to_kmer_altkmer:
#			ss_kmers = ss_to_kmer_altkmer[ss_name]
#			if pb_name not in pb_to_kmers:
#				pb_to_kmers[pb_name]=[set(),set()] # two lists corresponding to h0 and h1 kmers
#
#			pb_to_kmers[pb_name][h].add(ss_kmers[0])
#			pb_to_kmers[pb_name][1-h].add(ss_kmers[1])
#			
#			# TODO: compute pb interval....

# alternative
print("computing pb kmer intervals...")

with open(snakemake.input["SS_PB_minimap"]) as minimap:
	for line in minimap:
		if line.strip() == "":
			break
		sp=line.split()
		# skip the header line
		if sp[0] == "PBreadNames":
			continue

		pb_name=sp[0]
		ss_name=sp[1]
		lib_name=sp[2]
		ss_len=int(sp[6])
		ss_start=int(sp[7])
		ss_end=int(sp[8])
		strand=sp[9]
		pb_start=int(sp[14])
		pb_end=int(sp[15])
		cigar=sp[18]
		pb_clust=sp[19]
		# skip processing this line if the ss read is not clustered
		if ss_name not in ss_to_clust:
			continue
		ss_clust=ss_to_clust[ss_name]
		map_dir= -1 if strand=="-" else 1

		if (lib_name, ss_clust) not in lib_clust_to_haplo:
			# the SS lib is not WC in the cluster
			continue

		h = lib_clust_to_haplo[(lib_name, ss_clust)]
		if ss_name in ss_to_kmer_altkmer: # ss read covers a heterozygous position
			ss_kmers = ss_to_kmer_altkmer[ss_name]
			ss_het_pos=ss_to_kmer_pos_and_interval[ss_name][0]
			ss_kmer_interval=ss_to_kmer_pos_and_interval[ss_name][1]
			aln_dir = map_dir*ss_bubble_map_dir[ss_name]
			
			if ss_het_pos < ss_start or ss_het_pos > ss_end:
				# the aligned part of the ss read does not contain the het position
				continue


			# update ss_kmer_interval if it is not fully alligned to PB read
			ss_kmer_interval[0] = max(ss_kmer_interval[0], ss_start)
			ss_kmer_interval[1] = min(ss_kmer_interval[1], ss_end)

			if pb_name=='m160220_032405_00116_c100962052550000001823220307011691_s1_p0/101239/0_18961':
				print('ss_kmer_interval=', ss_kmer_interval)
			
			if aln_dir != 1: # either ss to bubble or ss to pb alignment direction is reverse
				ss_kmer_interval[0]=ss_len-ss_kmer_interval[0]-1
				ss_kmer_interval[1]=ss_len-ss_kmer_interval[1]-1
				ss_kmer_interval=(ss_kmer_interval[1], ss_kmer_interval[0])

			if pb_name=='m160220_032405_00116_c100962052550000001823220307011691_s1_p0/101239/0_18961':
				print('rc_ss_kmer_interval=', ss_kmer_interval)

			if pb_name not in pb_to_kmers:
				pb_to_kmers[pb_name]=[[],[]] # two lists corresponding to h0 and h1 kmers

			# computing the kmer interval in the pb read
			print('mapping index for pb read', pb_name, 'and ss read', ss_name)
			pb_kmer_interval_start = map_index(ss_kmer_interval[0]-ss_start, cigar)+pb_start
			pb_kmer_interval_end = map_index(ss_kmer_interval[1]-ss_start, cigar)+pb_start
			
			if pb_name=='m160220_032405_00116_c100962052550000001823220307011691_s1_p0/101239/0_18961':
				print('ss_kmer_interval=', ss_kmer_interval)
				print('ss_start=', ss_start)
				print('ss_kmer_interval-ss_start=', ss_kmer_interval[0]-ss_start, ss_kmer_interval[1]-ss_start)
				print('pb_start=', pb_start)
			
			if pb_kmer_interval_start=="" or pb_kmer_interval_start=="":
				# there are characters other than {M, X, I, D} in the cigar string
				#print("cigar = ", cigar, "there are characters other than {M, X, I, D} in the cigar string")
				continue

			pb_kmer_interval = (pb_kmer_interval_start, pb_kmer_interval_end)

			if aln_dir==1:
				pb_to_kmers[pb_name][h].append((ss_kmers[0], pb_kmer_interval))
				pb_to_kmers[pb_name][1-h].append((ss_kmers[1], pb_kmer_interval))
			else:
				pb_to_kmers[pb_name][h].append((reversecomp(ss_kmers[0]), pb_kmer_interval))
				pb_to_kmers[pb_name][1-h].append((reversecomp(ss_kmers[1]), pb_kmer_interval))

print("pb_to_kmers['m160220_032405_00116_c100962052550000001823220307011691_s1_p0/101239/0_18961']=", pb_to_kmers['m160220_032405_00116_c100962052550000001823220307011691_s1_p0/101239/0_18961'])
				

# reading pb fasta file and getting het kmer subsequences
print("outputting pb kmers...")

with open(snakemake.input["PB_fasta"]) as f:
	with open(snakemake.output[0], 'w') as out:
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
					for j in range(len(pb_to_kmers[name][0])):
						pb_kmer_interval = pb_to_kmers[name][0][j][1]
						subseq=seq[pb_kmer_interval[0]:(pb_kmer_interval[1]+1)]
						#print("name = ", name)
						#print("subseq = ", subseq)
						#print("h1_kmer = ", pb_to_kmers[name][0][j][0])
						#print("h2_kmer = ", pb_to_kmers[name][1][j][0])
						print(name + "\t" + subseq + "\t" + pb_to_kmers[name][0][j][0] + "\t" + pb_to_kmers[name][1][j][0], file=out)



#with open(snakemake.output[0], 'w') as out:
#	for pb in pb_to_kmers:
#		pb_kmers = pb_to_kmers[pb]
#		print(">" + pb + "\t#hetpos=" + str(len(pb_kmers[0])), file=out)
#		for h in range(2):
#			print(">haplotype" + str(h), file=out)
#			for kmer in pb_kmers[h]:
#				print(kmer, file=out)


