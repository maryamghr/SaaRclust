import pdb

from parsing import *


def phase_bubbles(ss_bubble_map_files, bubbles, ss_to_clust, lib_clust_to_haplo, bubble_phase_file, min_h0_frac=0.25, max_h0_frac=0.75):
	'''
	Computes the bubble allele corresponding to the first haplotype (h0) and writes it in the output file
	
	Args:
		ss_bubble_map_files: a list of name of files containing the exact mapping information of SS reads to SNV bubbles
		lib_clust_to_haplo:  A dictionaty that maps a pair (ss_lib, clust) to a haplotype (only for wc ss libs)		
		bubble_phase_file: The output file containing two columns for the bubble name and the bubble allele corresponding to the first haplotype, respectively
	'''
	
	#bubble_haplo_allele_count = {}
	#bubbles_with_invalid_alleles = {}
	#bubble_allele_to_unitig_name = {}
	
	unitig_haplo_cov = {}
	# for testing
	bubble_unitig_haplo_cov = {}

	if type(ss_bubble_map_files)==str:
		ss_bubble_map_files=[ss_bubble_map_files]

	for input_file in ss_bubble_map_files:
		with open(input_file) as f:
			for line in f:
				if line=="":
					break
				
				# skip the header line
				if line.startswith("SSname"):
					continue
					
				sp = line.split()
				
				ss_name, ss_lib, unitig_name, bubble_id, bubble_allele_id, is_reverse = sp[0], sp[1], sp[2], sp[3], sp[4], sp[5]

				if ss_name not in ss_to_clust:
					continue

				ss_clust = ss_to_clust[ss_name]

				if (ss_lib, ss_clust) not in lib_clust_to_haplo:
					continue

				haplo = lib_clust_to_haplo[(ss_lib, ss_clust)]
				
				#if unitig_name=="utg000518l":
				#	pdb.set_trace()
				
				if bubble_id=="None":
					# it's a non-bubble unitig
				#	if unitig_name=="utg000518l":
				#		print('hello in if')
					if unitig_name not in unitig_haplo_cov:
						unitig_haplo_cov[unitig_name]=[0,0]
						
					unitig_haplo_cov[unitig_name][haplo] +=1
					continue
					
				#if unitig_name=="utg000518l":
				#	print('hello')
				#	pdb.set_trace()
					
				# just for testing:
				if unitig_name not in bubble_unitig_haplo_cov:
					bubble_unitig_haplo_cov[unitig_name]=[0,0]
						
					bubble_unitig_haplo_cov[unitig_name][haplo] +=1

				# check invalid alleles...
				#if ss_name == "ERR1295715.948":
				#	pdb.set_trace()
				
				print('line =', line)
				assert(bubble_allele_id=='0' or bubble_allele_id=='1'), 'error in bubble '+bubble_id+': bubble allele is not binary'
				
				bubble = bubbles[int(bubble_id)]
				bubble_allele = bubble.allele0 if bubble_allele_id=='0' else bubble.allele1
				
				bubble_allele.haplo_coverage[haplo] += 1


	with open(bubble_phase_file, 'w') as out:
		print("unitigName\tbubbleName\tbubbleAllele\thaplotype", file=out)
		# printing bubble phase info
		for bubble_id, bubble in bubbles.items():
			
			h0_al0 = bubble.allele0.haplo_coverage[0]+bubble.allele1.haplo_coverage[1]
			h0_al1 = bubble.allele0.haplo_coverage[1]+bubble.allele1.haplo_coverage[0]

			#if h0_al0 == h0_al1:
			#	continue
				
			if h0_al0+h0_al1==0 or min_h0_frac < h0_al0/(h0_al0+h0_al1) < max_h0_frac:
				continue

			#pdb.set_trace()

			al0_haplo = 'H1' if h0_al0 > h0_al1 else 'H2'
			al1_haplo = 'H1' if h0_al0 < h0_al1 else 'H2'

			print(bubble.allele0.unitig_name + "\t" + str(bubble.id) + "\t0\t" + al0_haplo, file=out)
			print(bubble.allele1.unitig_name + "\t" + str(bubble.id) + "\t1\t" + al1_haplo, file=out)
			
		# printing unitig phase info
		#pdb.set_trace()
		
		for unitig_name, haplo_cov in unitig_haplo_cov.items():
			
			#if unitig_name=="utg000144l" or untigi_name=="utg000005l" or unitig_name=="utg000429l" or unitig_name=="utg000116l":
			#	pdb.set_trace()
			
			if sum(haplo_cov)==0 or min_h0_frac < haplo_cov[0]/sum(haplo_cov) < max_h0_frac:
				continue
			
			haplo = 'H1' if haplo_cov[0] > haplo_cov[1] else 'H2'
			print(unitig_name + "\tNone\tNone\t" + haplo, file=out)

