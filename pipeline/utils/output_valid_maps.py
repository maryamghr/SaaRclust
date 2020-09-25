import time
#!/usr/bin/python
from parsing import *
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


# Mapping ss read names to their lengths
def output_valid_maps(ss_fasta, ss_clust_file, clust_pair, bubble_fasta, mummer_files, out_file, log_file, libs, input_type="bubble", unitig_to_bubble_allele=None):
	'''
	input_type \in {"bubble, unitig"}
	unitig_to_bubble_allele: should be provided id input_type = "unitig"
	'''
	with open(log_file, 'w') as log:
		start_time = time.time()
		print('reading input files ...', file=log)

		if type(mummer_files)==str:
			mummer_files = [mummer_files]

		print('reading input files ...')
		start_time = time.time()
		
		bubbles_to_len = get_seq_len(bubble_fasta)
		
		print('elapsed time:', time.time()-start_time, 's')
		print('reading ss fasta files ...')
		start_time = time.time()
		
		ss_to_len = get_seq_len(ss_fasta)
		
		print('elapsed time:', time.time()-start_time, 's')
		print('reading ss clust files ...')
		start_time = time.time()
		
		ss_clust = get_ss_clust(ss_clust_file)
		
		print('elapsed time:', time.time()-start_time, 's')
		#print('reading ss haplo files ...')
		#start_time = time.time()
		
		#lib_clust_to_haplo = read_strandphaser_strand_states(phased_strand_states_file)

		print('number of processed ss reads =', len(ss_to_len), file=log)
		print('elapsed time:', time.time()-start_time, 's', file=log)

		start_time = time.time()
		print('processing the map file and outputting the valid maps ...', file=log)

		ss_map_count = {}
		ss_map_print = {}
		test_ss_map_count = {}

		#with open(out_file, 'w') as valid_map:
		#	print("SSname\tSSlib\tunitig_name\tbubbleName\tbubbleAllele\tisReverseMapped\tbubbleStart\tSSstart\talnLen", file=valid_map)
		for i in range(len(mummer_files)):
			start_time = time.time()
			print('processing', mummer_files[i])
			
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
						if ss_name not in ss_map_count:
							ss_map_count[ss_name] = 1
						else:
							ss_map_count[ss_name] +=1
						
						unitig_start = int(sp[1])
						ss_start = int(sp[2])
						aln_len = int(sp[3])
						unitig_len = bubbles_to_len[sp[0]]
							
						ss_lib_name = libs[i]
						ss_len = ss_to_len[ss_name]
						ss_cluster = ss_clust[ss_name] if ss_name in ss_clust else "None"
						#ss_haplo = lib_clust_to_haplo[(ss_lib_name, ss_cluster)] if (ss_lib_name, ss_cluster) in lib_clust_to_haplo else "None"

						#if sp[0] =="utg000220l" and ss_lib_name=="HG00733.P0IIL031" and aln_len > 100:
						#	pdb.set_trace()
						#	print('line =', line)
						#	if is_valid(unitig_start, unitig_len, ss_start, ss_len, aln_len):
						#		test_ss_map_count[ss_name] = ss_map_count[ss_name]
						#		
						#	print('test_ss_map_count =', test_ss_map_count)
						
						if ss_cluster not in clust_pair:
							continue
					
						if not is_valid(unitig_start, unitig_len, ss_start, ss_len, aln_len):
							continue

						if input_type=="bubble":
							bubble_sp = sp[0].split("_")
							bubble_num = bubble_sp[1]
							allele_num = str(int(bubble_sp[3])-1)

						else:
							bubble_num, allele_num = 'None', 'None'
							if sp[0] in unitig_to_bubble_allele:
								# the unitig is a bubble
								bubble_num, allele_num = unitig_to_bubble_allele[sp[0]]

						#print(ss_name + "\t" + ss_lib_name + "\t" + sp[0] + "\t" + bubble_num + "\t" + allele_num + "\t" + str(rev) + "\t" + str(unitig_start) + "\t" + str(ss_start) + "\t" + str(aln_len), file=valid_map)
						# chaning to deal with mummer bug (mummer sometimes outputs ss reads more than once)
						
						ss_map_print[ss_name] = ss_name + "\t" + ss_lib_name + "\t" + sp[0] + "\t" + bubble_num + "\t" + allele_num + "\t" + str(rev) + "\t" + str(unitig_start) + "\t" + str(ss_start) + "\t" + str(aln_len)
				
				#pdb.set_trace()
						
				print(time.time()-start_time, 's')
				start_time = time.time()

		print('outputing valid maps')
		with open(out_file, 'w') as valid_map:
			print("SSname\tSSlib\tunitig_name\tbubbleName\tbubbleAllele\tisReverseMapped\tbubbleStart\tSSstart\talnLen", file=valid_map) # remove ss_haplo later
			print("SSname\tSSlib\tbubbleName\tbubbleAllele\tisReverseMapped\tbubbleStart\tSSstart\talnLen", file=valid_map)
			for ss_name in ss_map_print:
				if ss_map_count[ss_name] > 1:
					continue
					
				print(ss_map_print[ss_name], file=valid_map)
				
		print(time.time()-start_time, 's')
			