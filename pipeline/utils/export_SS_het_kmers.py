import sys
import copy
#import zipfile
import gzip
#import re

# Note: q should be lower than the kmer length...
q = snakemake.params["het_kmer_len"]

print("q = ", q)

revcomp = {'a':'t', 'c':'g', 'g':'c', 't':'a', 'A':'T', 'C':'G', 'G':'C', 'T':'A'}

# returns the reverse complement of string seq
def reversecomp(seq):
	rc = ''
	for i in range(len(seq)):
		rc = revcomp[seq[i]] + rc
	return rc



# This function adds the set of heterozygous kmers for the two chains of a bubble
# bubble_allele_to_kmers is a dictionary that maps a bubble chain to the set of its heterozygous kmers
# seq is a list of two strings corresponding to the two chains of the bubble with the id bubble_id
# q is the length of the heterozygous kmer to be added to the dictionary

def add_bubble_kmers(bubble_allele_to_kmers, bubble_id, seq, q):
	for al in range(2):
		bubble_allele_to_kmers[(bubble_id, al)]={}
	d = []
	for i in range(len(seq[0])):
		if seq[0][i] != seq[1][i]:
			if i < q:
				print('bubble', bubble_id)
				print(seq[0])
				print(seq[1])
				sys.exit('Error in bubble: the heterozygous position is less than q!!!')
			for al in range(2):
				bubble_allele_to_kmers[(bubble_id, al)][i]=seq[al][i-q:(i+q+1)]


# This function returns the index of a string which is aligned to the idx position of another string
def map_index(idx, cigar):
	# counting the number of insertions and deletions that happen before each index in the aligned text	
	num_insertion=[]
	num_deletion= []
	insertion = 0
	deletion = 0
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
			if allele > 1:
				bubbles_with_invalid_alleles.add(bubble_id)

			if bubble_id != prev_bubble_id:			
				if prev_clust!="None" and num_bubble_rep == 2 and prev_bubble_id not in bubbles_with_invalid_alleles:				
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
			num_bubble_rep = num_bubble_rep + 1
			if allele < 2:				
				seq[allele] = line.strip()

	add_bubble_kmers(bubble_allele_to_kmers, bubble_id, seq, q)


print("bubble_allele_to_kmers:")
print(bubble_allele_to_kmers)


ss_to_kmer_altkmer={}
ss_to_kmer_pos_and_interval={}
ss_to_clust={}
ss_bubble_map_dir={}


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
				
				if (bubble_id, allele) not in bubble_allele_to_kmers:
					continue

				# checking which het pos in the bubble is covered by the ss read
				bubble_start = int(sp[8])-1
				ss_start = int(sp[9])-1
				aln_len = int(sp[10])
				ss_end = ss_start + aln_len -1 # 0-based

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
						if ss_kmer_interval[1] > ss_end:	 #aln_len + ss_start-2:
							# het pos is at the end of the ss read
							# cut the beginning of kmer and alt kmer
							cut = ss_kmer_interval[1] - ss_end 	#- aln_len + ss_start
							kmer = kmer[:-cut]
							alt_kmer = alt_kmer[:-cut]
							ss_kmer_interval[1] = ss_end	 #aln_len + ss_start-2
							#ss_kmer_interval = [ss_kmer_interval[0], aln_len + ss_start-1]

						ss_to_kmer_altkmer[ss_name]=(kmer, alt_kmer)
						ss_to_kmer_pos_and_interval[ss_name]=[ss_het_pos, ss_kmer_interval]
						ss_to_clust[ss_name]=ss_clust
						ss_bubble_map_dir[ss_name]=map_dir
						break


print("ss_to_kmer_altkmer =", ss_to_kmer_altkmer)
print("ss_to_kmer_pos_and_interval =", ss_to_kmer_pos_and_interval)
print("ss_to_clust =", ss_to_clust)
print("ss_bubble_map_di r=", ss_bubble_map_dir)





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


print("lib_clust_to_haplo =", lib_clust_to_haplo)


# pb_to_kmers is a dictionary that maps each pb read name to a list of two lists: the first list contains h0 kmers with their intervals and the second one contains h1 kmers with their intervals in the pb read
pb_to_kmers={}


print("computing pb kmer intervals...")

#with zipfile.ZipFile(snakemake.input["SS_PB_minimap"]) as z:
with gzip.open(snakemake.input["SS_PB_minimap"], 'rb') as minimap:
	for line in minimap:
		line = line.decode("utf-8")
		if line.strip() == "":
			break
		sp=line.split()
		# skip the header line
		if sp[0] == "PBreadNames":
			continue

		#print(line.decode("utf-8"), 'class =', type(line.decode("utf-8")))
		pb_name=sp[0]
		ss_name=sp[1]
		lib_name=sp[2]
		ss_len=int(sp[6])
		ss_start=int(sp[7])
		ss_end=int(sp[8])
		strand=sp[9]
		pb_start=int(sp[14])
		pb_end=int(sp[15])
		aln_len=int(sp[17])
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
			ss_kmers = copy.copy(ss_to_kmer_altkmer[ss_name])
			ss_het_pos = copy.copy(ss_to_kmer_pos_and_interval[ss_name][0])
			ss_kmer_interval = copy.copy(ss_to_kmer_pos_and_interval[ss_name][1])
			aln_dir = map_dir*copy.copy(ss_bubble_map_dir[ss_name])

		
			if aln_dir != 1: # either ss to bubble or ss to pb alignment direction is reverse
				ss_het_pos=ss_len-ss_het_pos-1
				ss_kmer_interval[0] = ss_len - ss_kmer_interval[0]-1
				ss_kmer_interval[1] = ss_len - ss_kmer_interval[1]-1
				ss_kmer_interval.reverse()

			if ss_het_pos < ss_start or ss_het_pos > ss_end-1:
				# the aligned part of the ss read does not contain the het position
				continue


			# update ss_kmer_interval if it is not fully alligned to PB read
			ss_kmer_interval[0] = max(ss_kmer_interval[0], ss_start)
			ss_kmer_interval[1] = min(ss_kmer_interval[1], ss_end-1)

			if pb_name not in pb_to_kmers:
				pb_to_kmers[pb_name]=[[],[]] # two lists corresponding to h0 and h1 kmers

			# computing the kmer interval in the pb read
		
			pb_kmer_interval_start = map_index(ss_kmer_interval[0]-ss_start, cigar)+pb_start
			pb_kmer_interval_end = map_index(ss_kmer_interval[1]-ss_start, cigar)+pb_start
		
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


print("pb_to_kmers =", pb_to_kmers)


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
						print(name + "\t" + subseq + "\t" + pb_to_kmers[name][0][j][0] + "\t" + pb_to_kmers[name][1][j][0], file=out)

