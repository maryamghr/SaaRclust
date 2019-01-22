q = 10

bubble_allele_to_kmers = {}

with open(snakemake.input["bubbles"]) as bubbles:
	prev_bubble_id=-1
	seq=["", ""]
	for line in bubbles:
		if line!="" and line[0]==">":
			bubble_name=line.split()[0][1:]
			sp=bubble_name.split("_")
			bubble_id=int(sp[1])
			allele=int(sp[3])-1
			if bubble_id != prev_bubble_id:
				# process the previous bubble
				if len(seq[0])!=len(seq[1]):
					print("warning: the sequence lengths are different!!!")
					break
				bubble_allele_to_kmers[(bubble_id, allele)]={}
				d = []
				for i in range(len(seq[0])):
					if seq[0][i] != seq[1][i]:
						for al in range(2):
							bubble_allele_to_kmers[(bubble_id, al)][i]=seq[al][i-q:(i+q+1)]
						
				
				prev_bubble_id = bubble_id
				seq=["", ""]
			else:
				seq[allele] = line.strip()

with open(snakemake.input["SS_bubble_map"]) as f:
	ss_to_kmer_altkmer={}
	for line in f:
		if line!="" and line[0]==">":
			ss_name=line.split()[1]
		else:
			sp=line.split()
			bubble_id=int(sp[0].split("_")[1])
			allele=int(sp[0].split("_")[3])
			
			if (bubble_id, allele) not in bubble_allele_to_kmers:
				continue
			# checking which het pos in the bubble in covered by the ss read
			bubble_start = int(sp[1])-1
			aln_len = int(sp[3])
			for het_pos in bubble_allele_to_kmers:
				if het_pos >= bubble_start and het_pos - bubble_start < aln_len:
					# heterozygous position het_pos is covered by the strand-seq read
					ss_to_kmer_altkmer[ss_name]=(bubble_allele_to_kmers[(bubble_id, allele)][het_pos], bubble_allele_to_kmers[(bubble_id, 1-allele)][het_pos])
					break

with open(snakemake.input["SS_haplo_strand_state"]) as states:
	ss_state[cluster_name]=0/1

h = ss_state[snakemake.wildcards["cluster"]]

with open(snakemake.input["SS_PB_minimap"]) as minimap:
	pb_to_kmers={}
	for lines in minimap:
		sp=lines.split()
		pb_name=sp[0]
		ss_name=sp[1]
		if ss_name in ss_to_kmer_altkmer:
			ss_kmers = ss_to_kmer_altkmer[ss_name]
			if pb_name not in pb_to_kmers:
				pb_to_kmers[pb_name]=[set(),set()] # two lists corresponding to h0 and h1 kmers

			pb_to_kmers[pb_name][h].add(ss_kmers[0])
			pb_to_kmers[pb_name][1-h].add(ss_kmers[1])

with open(snakemake.output[0], 'w') as out:
	for pb in pb_to_kmers:
		print(">", pb)
		for h in range(2):
			for 
