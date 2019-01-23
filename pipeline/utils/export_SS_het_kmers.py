import sys

q = snakemake.params["het_kmer_len"]

def add_bubble_kmers(bubble_allele_to_kmers, bubble_id, seq):
	for al in range(2):
		bubble_allele_to_kmers[(prev_bubble_id, al)]={}
	d = []
	for i in range(len(seq[0])):
		if seq[0][i] != seq[1][i]:
			for al in range(2):
				bubble_allele_to_kmers[(prev_bubble_id, al)][i]=seq[al][i-q:(i+q+1)]

bubble_allele_to_kmers = {}

with open(snakemake.input["bubbles"]) as bubbles:
	prev_bubble_id=-1
	seq=["", ""]
	for line in bubbles:
		if line=="":
			break

		if line[0]==">":
			bubble_name=line.split()[0][1:]
			sp=bubble_name.split("_")
			bubble_id=int(sp[1])
			allele=int(sp[3])-1
			if bubble_id != prev_bubble_id:
				# process the previous bubble
				if len(seq[0])!=len(seq[1]):
					print(seq)
					sys.exit("Error: the sequence lengths are different!!!")
				
				if prev_bubble_id != -1:
					add_bubble_kmers(bubble_allele_to_kmers, prev_bubble_id, seq)		
				
				prev_bubble_id = bubble_id
				seq=["", ""]
		else:
			seq[allele] = line.strip()

	add_bubble_kmers(bubble_allele_to_kmers, prev_bubble_id, seq)

### test # passed
print("bubble_allele_to_kmers:")
print(bubble_allele_to_kmers)


ss_to_kmer_altkmer={}

with open(snakemake.input["SS_bubble_map"]) as f:
	for line in f:
		if line=="":
			break
		
		if line[0]==">":
			ss_name=line.split()[1].split("_lib")[0]
		else:
			sp=line.split()
			bubble_id=int(sp[0].split("_")[1])
			allele=int(sp[0].split("_")[3])-1
			print("ss_name = ", ss_name, "bubble_id = ", bubble_id, ", allele = ", allele)
			
			if (bubble_id, allele) not in bubble_allele_to_kmers:
				continue
			# checking which het pos in the bubble in covered by the ss read
			bubble_start = int(sp[1])-1
			aln_len = int(sp[3])
			for het_pos in bubble_allele_to_kmers[(bubble_id, allele)]:
				if het_pos >= bubble_start and het_pos - bubble_start < aln_len:
					# heterozygous position het_pos is covered by the strand-seq read
					ss_to_kmer_altkmer[ss_name]=(bubble_allele_to_kmers[(bubble_id, allele)][het_pos], bubble_allele_to_kmers[(bubble_id, 1-allele)][het_pos])
					break


### test
print("ss_to_kmer_altkmer:")
print(ss_to_kmer_altkmer)

# mapping (lib, clust) pairs to haplotype
lib_clust_to_haplo = {}

for f in snakemake.input["SS_haplo_strand_states"]:
	file_name = f.split("/")[-1]
	lib_name=file_name.split("_haplo_strand_states.data")[0]
	with open(f) as states:
		for line in states:
			# process the line if it is not empty nor the header line
			if line!="" and line[0]=="V":
				sp = line.split()
				lib_clust_to_haplo[(lib_name, sp[0])]=int(sp[1])-1

### test
#print("lib_clust_to_haplo:")
#print(lib_clust_to_haplo)


cluster = snakemake.wildcards["cluster"]
print("cluster = ", cluster)

pb_to_kmers={}

with open(snakemake.input["SS_PB_minimap"]) as minimap:
	for line in minimap:
		if line.strip() == "":
			break
		sp=line.split()
		pb_name=sp[0]
		ss_name=sp[1]
		lib_name=sp[2]
		if (lib_name, cluster) not in lib_clust_to_haplo:
			# the SS lib is not WC in the cluster
			continue
		h = lib_clust_to_haplo[(lib_name, cluster)]
		if ss_name in ss_to_kmer_altkmer:
			ss_kmers = ss_to_kmer_altkmer[ss_name]
			if pb_name not in pb_to_kmers:
				pb_to_kmers[pb_name]=[set(),set()] # two lists corresponding to h0 and h1 kmers

			pb_to_kmers[pb_name][h].add(ss_kmers[0])
			pb_to_kmers[pb_name][1-h].add(ss_kmers[1])


with open(snakemake.output[0], 'w') as out:
	for pb in pb_to_kmers:
		pb_kmers = pb_to_kmers[pb]
		print(">" + pb + "\t#hetpos=" + str(len(pb_kmers[0])), file=out)
		for h in range(2):
			print(">haplotype" + str(h), file=out)
			for kmer in pb_kmers[h]:
				print(kmer, file=out)


