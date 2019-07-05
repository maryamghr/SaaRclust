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


with open(snakemake.output[0], 'w') as valid_map:
	with open(snakemake.input[0]) as f:
		print("SSname\tSSlib\tSSclust\tbubbleName\tbubbleAllele\tbubbleChrom\tbubbleFlag\tisReverseMapped\tbubbleStart\tSSstart\talnLen", file=valid_map)
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
				unitig_len = int(sp[0].split("_len:")[-1])
				ss_name_sp = ss_name.split("_")

				if (not ss_name_sp[4].isdigit()) or ss_name_sp[5] == "clust:None":
					continue # the ground true chromosome is unknown (corresponding to the chrom names which contain "_"), or the ss read is not clustered
				
				ss_read_name = ss_name_sp[0]
				ss_lib_name = ss_name_sp[1] + "_" + ss_name_sp[2]
				ss_clust = ss_name_sp[5].split("clust:")[1]
				ss_len = int(ss_name_sp[6].split("len:")[1])

				if is_valid(unitig_start, unitig_len, ss_start, ss_len, aln_len):
					bubble_sp = sp[0].split("_")
					bubble_num = bubble_sp[1]
					allele_num = str(int(bubble_sp[3])-1)
					bubble_flag = bubble_sp[6]
					bubble_chrom = bubble_sp[7]
					print(ss_read_name + "\t" + ss_lib_name + "\t" + ss_clust + "\t" + bubble_num + "\t" + allele_num + "\t" + bubble_chrom + "\t" + bubble_flag + "\t" + str(rev) + "\t" + str(unitig_start) + "\t" + str(ss_start) + "\t" + str(aln_len), file=valid_map)
