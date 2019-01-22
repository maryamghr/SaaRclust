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

with open(snakemake.output["matrix"], 'w') as out:
	with open(snakemake.output["map"], 'w') as valid_map:
		with open(snakemake.input[0]) as f:
			for line in f:
				sp = line.split()
				if len(sp) == 0:
					break
				if sp[0] == ">":
					ss_name =  sp[1]
				else: # reading snv bubble information: there is an exact match with a bubble
					unitig_start = int(sp[1])
					ss_start = int(sp[2])
					aln_len = int(sp[3])
					unitig_len = int(sp[0].split("_len:")[-1])
					ss_name_sp = ss_name.split("_len:")
					ss_len = int(ss_name_sp[-1])
					libname = ss_name_sp[0].split("_lib:")[-1]
					libs.add(libname)

					if is_valid(unitig_start, unitig_len, ss_start, ss_len, aln_len):
						print("> " + ss_name + "\n" + line.strip(), file=valid_map)
						bubble_sp = sp[0].split("_")
						bubble_num = bubble_sp[1]
						allele_num = str(int(bubble_sp[3])-1)
						if bubble_num in bubble_cov:
							bubble_cov[bubble_num][libname] = allele_num
						else:
							bubble_cov[bubble_num] = {libname : allele_num}
	# mapping libnames to their indices
	libs = list(libs)
	lib_name_to_idx = {}
	for i in range(len(libs)):
		lib_name_to_idx[libs[i]] = i
		
	# writing the counts in the output
	# print header
	print("bubble\\cell\t" + "\t".join(libs), file=out)

	# print bubble coverages
	gap = ['-']*len(libs)
	for bubble in bubble_cov:
		b = gap
		for cell in bubble_cov[bubble]:
			b[lib_name_to_idx[cell]] = bubble_cov[bubble][cell]
		print(bubble + "\t" + "\t".join(b), file=out)
