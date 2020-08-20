import time
#!/usr/bin/python
import pdb


def is_valid(s1, l1, s2, l2, aln_len):
	if aln_len == min(l1, l2): # one string contains the other
		return (True)
	elif s1+aln_len-1 == l1 and s2 == 1: # suffix of str1 matches prefix of str2
		return (True)
	elif s2+aln_len-1 == l2 and s1 == 1: # suffix of str2 matches prefix of str1
		return (True)

	return (False)


def print_dict_head(dic, num):
	for i in range(num):
		key = list(dic.keys())[i]
		print(key, ':', dic[key])


def get_seq_len(fasta_files):
	if type(fasta_files)==str:
		fasta_files = [fasta_files]

	seq_to_len = {}
	for fasta_file in fasta_files:
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


def map_unitig_to_bubble_allele(bubble_file):
	unitig_to_bubble_allele = {}
	with open(bubble_file) as bubbles:
		for line in bubbles:
			if line == "":
				break

			sp = line.strip().split("_")
			bubble_id, bubble_allele, unitig = sp[1], str(int(sp[3])-1), sp[-1]

			unitig_to_bubble_allele[unitig] = (bubble_id, bubble_allele)

	return unitig_to_bubble_allele


# Mapping ss read names to their lengths
def output_valid_maps(ss_fasta, bubble_fasta, mummer_files, out_file, log_file, libs, input_type="bubble", unitig_to_bubble_allele=None):
	'''
	input_type \in {"bubble, unitig"}
	unitig_to_bubble_allele: should be provided id input_type = "unitig"
	'''
	with open(log_file, 'w') as log:
		start_time = time.time()
		print('reading the ss fastq file and mapping names to len ...', file=log)

		if type(mummer_files)==str:
			mummer_files = [mummer_files]

		bubbles_to_len = get_seq_len(bubble_fasta)
		ss_to_len = get_seq_len(ss_fasta)

		print('number of processed ss reads =', len(ss_to_len), file=log)
		print('elapsed time:', time.time()-start_time, 's', file=log)

		start_time = time.time()
		print('processing the map file and outputting the valid maps ...', file=log)

		with open(out_file, 'w') as valid_map:
			print("SSname\tSSlib\tbubbleName\tbubbleAllele\tisReverseMapped\tbubbleStart\tSSstart\talnLen", file=valid_map)
			for i in range(len(mummer_files)):
				with open(mummer_files[i]) as f:
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

							ss_lib_name = libs[i]
							ss_len = ss_to_len[ss_name]

							if not is_valid(unitig_start, unitig_len, ss_start, ss_len, aln_len):
								continue

							if input_type=="bubble":
								bubble_sp = sp[0].split("_")
								bubble_num = bubble_sp[1]
								allele_num = str(int(bubble_sp[3])-1)

							else:
								bubble_num, allele_num = sp[0], 'None'
								if sp[0] in unitig_to_bubble_allele:
									# the unitig is a bubble
									bubble_num, allele_num = unitig_to_bubble_allele[sp[0]]

							
							print(ss_read_name + "\t" + ss_lib_name + "\t" + bubble_num + "\t" + allele_num + "\t" + str(rev) + "\t" + str(unitig_start) + "\t" + str(ss_start) + "\t" + str(aln_len), file=valid_map)
						
