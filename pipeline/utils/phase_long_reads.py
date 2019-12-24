def phase_long_reads(bubble_pb_kmer_file, pb_phase_file):
	pb_to_h0_h1_dist = {}
	pb_name_to_num_kmers = {}
	pb_name_to_num_phased_kmers = {}

	with open(bubble_pb_kmer_file) as f:
		# skip the header line
		next(f)
		for line in f:
			if line=="":
				break

			sp = line.split()
			bubble_id, bubble_allele, pb_name, haplo, kmer_edit_dist = int(sp[0]), int(sp[1]), sp[2], sp[5], int(sp[6])
			if pb_name not in pb_name_to_num_kmers:
				pb_name_to_num_kmers[pb_name] = 0
				pb_name_to_num_phased_kmers[pb_name] = 0

			pb_name_to_num_kmers[pb_name] += 1

			if pb_name not in pb_to_h0_h1_dist:
				pb_to_h0_h1_dist[pb_name] = [0,0]

			if haplo == "?":
				# the bubble is not phased
				continue

			pb_name_to_num_phased_kmers[pb_name] += 1

			h = int(haplo)

			pb_to_h0_h1_dist[pb_name][h] += kmer_edit_dist

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
