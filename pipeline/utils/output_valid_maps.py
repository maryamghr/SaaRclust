import time
#!/usr/bin/python


def is_valid(s1, l1, s2, l2, aln_len):
	if aln_len == min(l1, l2): # one string contains the other
		return (True)
	elif s1+aln_len-1 == l1 and s2 == 1: # suffix of str1 matches prefix of str2
		return (True)
	elif s2+aln_len-1 == l2 and s1 == 1: # suffix of str2 matches prefix of str1
		return (True)

	return (False)

ss_name = ""
bubble_cov = {}
libs = set()

def print_dict_head(dic, num):
	for i in range(num):
		key = list(dic.keys())[i]
		print(key, ':', dic[key])


def get_seq_len(fasta_file):
	seq_to_len = {}
	with open(fasta_file) as f:
		seq_name = ""
		for line in f:
			if line == "":
				break

			if line.startswith('>'):
				seq_name = line.strip()[1:]

			else:
				seq_to_len[seq_name] = len(line.strip())

	return (seq_to_len)


# Mapping ss read names to their lengths
with open(snakemake.log[0], 'w') as log:
	start_time = time.time()
	print('reading the ss fastq file and mapping names to len ...', file=log)

	ss_to_len = get_seq_len(snakemake.input['ss_reads'])
	bubbles_to_len = get_seq_len(snakemake.input['bubbles'])

	print('number of processed ss reads =', len(ss_to_len), file=log)
	print('elapsed time:', time.time()-start_time, 's', file=log)

	start_time = time.time()
	print('processing the map file and outputting the valid maps ...', file=log)

	with open(snakemake.output[0], 'w') as valid_map:
		with open(snakemake.input['map']) as f:
			print("SSname\tSSlib\tbubbleName\tbubbleAllele\tisReverseMapped\tbubbleStart\tSSstart\talnLen", file=valid_map)
			rev = False
			for line in f:
				sp = line.split()
				if len(sp) == 0:
					break
				if sp[0] == ">":
					ss_name =  sp[1]
					if len(sp)==3 and sp[2]=="Reverse":
						rev=True
					else:
						rev=False

				else: # reading snv bubble information: there is an exact match with a bubble
					unitig_start = int(sp[1])
					ss_start = int(sp[2])
					aln_len = int(sp[3])
					unitig_len = bubbles_to_len[sp[0]]
					ss_read_name = ss_name

					ss_lib_name = snakemake.wildcards['x']
					ss_len = ss_to_len[ss_name]

					if is_valid(unitig_start, unitig_len, ss_start, ss_len, aln_len):
						bubble_sp = sp[0].split("_")
						bubble_num = bubble_sp[1]
						allele_num = str(int(bubble_sp[3])-1)
						print(ss_read_name + "\t" + ss_lib_name + "\t" + bubble_num + "\t" + allele_num + "\t" + str(rev) + "\t" + str(unitig_start) + "\t" + str(ss_start) + "\t" + str(aln_len), file=valid_map)
						
