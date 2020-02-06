import time

def phase_long_reads(bubble_pb_kmer_file, bubble_phase_file, pb_phase_file):
	bubble_to_allele0_haplo = {}
	pb_to_h0_h1_dist = {}
	pb_name_to_num_kmers = {}
	pb_name_to_num_phased_kmers = {}

	with open(bubble_phase_file) as bubble_phase:
		# skip header
		next(bubble_phase)
		for line in bubble_phase:
			if line=="":
				break

			sp = line.split()
			bubble_to_allele0_haplo[sp[0]] = int(sp[1])

	with open(bubble_pb_kmer_file) as f:
		# skip the header line
		next(f)
		for line in f:
			if line=="":
				break

			sp = line.split()
			bubble_id, bubble_allele, pb_name, kmer_edit_dist = sp[0], int(sp[1]), sp[2], int(sp[5])
			if pb_name not in pb_name_to_num_kmers:
				pb_name_to_num_kmers[pb_name] = 0
				pb_name_to_num_phased_kmers[pb_name] = 0

			pb_name_to_num_kmers[pb_name] += 1

			if pb_name not in pb_to_h0_h1_dist:
				pb_to_h0_h1_dist[pb_name] = [0,0]

			if bubble_id not in bubble_to_allele0_haplo:
				# the bubble is not phased
				continue

			pb_name_to_num_phased_kmers[pb_name] += 1

			allele0_haplo = bubble_to_allele0_haplo[bubble_id]

			haplo = allele0_haplo if bubble_allele==0 else 1-allele0_haplo

			pb_to_h0_h1_dist[pb_name][haplo] += kmer_edit_dist

	with open(pb_phase_file, 'w') as out:
		print("PBname\tnum_het_kmers\tnum_phased_het_mers\th0_edit_dist\th1_edit_dist\thaplotype", file=out)
		for pb_name in pb_to_h0_h1_dist:
			haplo_edit_dist = pb_to_h0_h1_dist[pb_name]
			haplotype = "?"
			
			if haplo_edit_dist[0] < haplo_edit_dist[1]:
				haplotype = "0"
			elif haplo_edit_dist[1] < haplo_edit_dist[0]:
				haplotype = "1"
		
			print(pb_name + "\t" + str(pb_name_to_num_kmers[pb_name]) + "\t" + str(pb_name_to_num_phased_kmers[pb_name]) + "\t" + str(haplo_edit_dist[0]) + "\t" + str(haplo_edit_dist[1]) + "\t" + haplotype, file=out)



def update_bubble_phase(bubble_pb_kmer_file_list, pb_phase_file_list, bubble_phase_file):
	pb_to_haplo = {}
	# bubble_alleles_to_haplotypes_edit_dist maps each bubble to a four-tuple: (d(allele0, haplo0),  d(allele0, haplo1),  d(allele1, haplo0),  d(allele1, haplo1))
	bubble_alleles_to_haplotypes_edit_dist = {}
	
	for pb_phase_file in pb_phase_file_list:
		with open(pb_phase_file) as phase_file:
			# skip the header line
			next(phase_file)
			for line in phase_file:
				if line=="":
					break

				sp = line.split()

				if sp[-1] != '?':
					pb_to_haplo[sp[0]] = sp[-1]

	
	for bubble_pb_kmer_file in bubble_pb_kmer_file_list:
		print('processing file', bubble_pb_kmer_file)
		start_time = time.time()
		with open(bubble_pb_kmer_file) as kmer_file:
			# skip the header line
			next(kmer_file)
			for line in kmer_file:
				if line=="":
					break

				sp = line.split()
				bubble_id, bubble_allele, pb_name, kmer_edit_dist = int(sp[0]), sp[1], sp[2], int(sp[5])
				if pb_name not in pb_to_haplo:
					continue

				pb_haplo = pb_to_haplo[pb_name]

				if bubble_id not in bubble_alleles_to_haplotypes_edit_dist:
					bubble_alleles_to_haplotypes_edit_dist[bubble_id] = [0,0,0,0]

				# finding the corresponding index in the allele_haplo_edit dist
				binary_index = bubble_allele + pb_haplo
				index = int(binary_index, 2)

				bubble_alleles_to_haplotypes_edit_dist[bubble_id][index] += kmer_edit_dist
		print('elapsed time =', time.time()-start_time, 's')

	with open(bubble_phase_file, 'w') as out:
		print("bubbleName\thaplotype0Allele", file=out)
		for bubble_id in bubble_alleles_to_haplotypes_edit_dist:
			dist0 = bubble_alleles_to_haplotypes_edit_dist[bubble_id][int('00',2)]+bubble_alleles_to_haplotypes_edit_dist[bubble_id][int('11',2)]
			dist1 = bubble_alleles_to_haplotypes_edit_dist[bubble_id][int('01',2)]+bubble_alleles_to_haplotypes_edit_dist[bubble_id][int('10',2)]

			if dist0 == dist1:
				continue
			
			bubble_haplo = 0 if dist0 < dist1 else 1
		
			print(bubble_id, "\t", bubble_haplo, file=out)


