#from __future__ import division
import gzip

min_overlap = snakemake.params[0]

with gzip.open(snakemake.input["minimap"]) as minimap:
	true = 0
	false = 0

	for line in minimap:
		line = line.decode()
		if line.strip() == "":
			break

		sp=line.split()
		# skip the header line
		if sp[0] == "SSreadNames":
			continue

		lib_name = sp[1]
		ss_flag = sp[2]
		ss_chrom = sp[3]
		ss_pos = int(sp[4])
		ss_len = int(sp[5])
		ss_start = int(sp[6])
		ss_end = int(sp[7])-1
		strand = sp[8]
		pb_flag = sp[10]
		pb_chrom = sp[11]
		pb_pos = int(sp[12])
		pb_len = int(sp[13])
		pb_start = int(sp[14])
		pb_end = int(sp[15])-1

		if ss_flag.isdigit():
			ss_flag = int(ss_flag)
		else:
			continue

		if pb_flag.isdigit():
			pb_flag = int(pb_flag)
		else:
			continue

		ss_end_pos = ss_pos + ss_len - 1
		pb_end_pos = pb_pos + pb_len - 1
		
		# checking whether the alignment is true or false
		# comparing chroms
		if ss_chrom != pb_chrom:
			false = false+1
			continue

		# computing SS and PB mapping directions
		ss_binary_flag = bin(ss_flag)[2:]
		pb_binary_flag = bin(pb_flag)[2:]

		# 1 stands for forward direction and -1 for backward direction
		ss_dir = 1
		pb_dir = 1
		
		if len(ss_binary_flag)>4 and ss_binary_flag[-5]=='1':
			ss_dir = -1

		if len(pb_binary_flag)>4 and pb_binary_flag[-5]=='1':
			pb_dir = -1

		aln_dir = 1 if strand=='+' else -1

		# checking whether the alignment direction is right
		if ss_dir * pb_dir * aln_dir == -1:
			false = false + 1
			continue

		# checking whether the strand-seq interval is completely inside PB interval in the reference genome
		if ss_pos > pb_pos and ss_end_pos < pb_end_pos:
			true = true + 1
		else:
			false = false + 1

		# computing the ss to pb aln interval in the reference genome
		#aln_start_in_ref = -1
		#aln_end_in_ref   = -1

		#if pb_dir == 1:
		#	aln_start_in_ref = pb_pos + pb_start
		#	aln_end_in_ref   = pb_pos + pb_end
		#else:
		#	aln_start_in_ref = pb_end_pos - pb_end
		#	aln_end_in_ref   = pb_end_pos - pb_start

		#mapped_ss_start_in_ref = -1
		#mapped_ss_end_in_ref   = -1

		#if ss_dir == 1:
		#	mapped_ss_start_in_ref = ss_pos + ss_start
		#	mapped_ss_end_in_ref   = ss_pos + ss_end
		#else:
		#	mapped_ss_start_in_ref = ss_end_pos - ss_end
		#	mapped_ss_end_in_ref   = ss_end_pos - ss_start

		# checking whether the two intervals have enough overlap
		#if mapped_ss_end_in_ref < aln_start_in_ref or aln_end_in_ref < mapped_ss_start_in_ref:
			# the two intervals are completely disjoint
		#	false = false+1
		#	continue
		
		#overlap_len = max(mapped_ss_end_in_ref, aln_end_in_ref) - max(mapped_ss_start_in_ref, aln_start_in_ref)
		#max_interval_len = max(mapped_ss_end_in_ref-mapped_ss_start_in_ref+1, aln_end_in_ref-aln_start_in_ref+1)

		#if overlap_len < max_interval_len*min_overlap:
			# the two intervals do not have enough overlap
		#	false = false+1
		#	continue

f = open(snakemake.input["log"], 'r').readlines()

user_time=f[-2].split()[1]
user_time=user_time.split('m')
sys_time =f[-1].split()[1]
sys_time=sys_time.split('m')

user_t = float(user_time[0])*60 + float(user_time[1][:-1])
sys_t  = float(sys_time[0])*60  + float(sys_time[1][:-1])

with open(snakemake.output[0], 'w') as out:
	print("num true alignments\t"  + str(true),  file=out)
	print("num false alignments\t" + str(false), file=out)
	#print("fraction of true alignments\t" + str(true/(true+false)), file=out)
	print("minimap alignment CPU time\t" + str(user_t + sys_t) + 's', file=out)
