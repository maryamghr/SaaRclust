from __future__ import division

def get_ss_counts_in_chroms(ss_counts_file):

	ss_counts_in_chroms = {}

	with open(ss_counts_file) as f:
		# skip the header
		next(f)

		for line in f:
			if line=="":
				break

			chrom, w, c, w_h0, w_h1, c_h0, c_h1 = line.split()

			ss_counts_in_chroms[chrom] = {'w':int(w), 'c':int(c), 'w_h0':int(w_h0), 'w_h1':int(w_h1), 'c_h0':int(c_h0), 'c_h1':int(c_h1)}

	return ss_counts_in_chroms

def get_chrom_dir_to_clust_map(clust_partners_file):
	chrom_dir_to_clust = {}

	with open(clust_partners_file) as f:
		# skip the header line
		next(f)
		for line in f:
			if line == "":
				break

			sp = line.split()
			chrom_dir_to_clust[sp[0]] = sp[1]

	return chrom_dir_to_clust



def output_ss_haplo_strand_states(chrom_dir_to_clust, ss_counts_in_chroms, ss_haplo_strand_states_file, min_w_frac_in_wc_state=0.4, max_w_frac_in_wc_state=0.6, max_haplo_count_ratio=0.15):

	with open(ss_haplo_strand_states_file, 'w') as out:
		print('cluster\thaplotype\tchromosome_direction', file=out)

		for chrom in ss_counts_in_chroms:
			ss_chrom_counts = ss_counts_in_chroms[chrom]
			w, c, w_h0, w_h1, c_h0, c_h1 = ss_chrom_counts['w'], ss_chrom_counts['c'], ss_chrom_counts['w_h0'], ss_chrom_counts['w_h1'], ss_chrom_counts['c_h0'], ss_chrom_counts['c_h1']

			if w_h0 + w_h1 + c_h0 + c_h1 == 0:
				# there are no haplotaged ss reads from this chrom
				continue

			w_frac = w/(w+c)

			if w_frac < min_w_frac_in_wc_state or w_frac > max_w_frac_in_wc_state:
				# the lib is not of type wc in this chrom
				continue

			w_h0_c_h1_count = w_h0 + c_h1
			w_h1_c_h0_count = w_h1 + c_h0

			haplo_count_ratio = min(w_h0_c_h1_count, w_h1_c_h0_count)/max(w_h0_c_h1_count, w_h1_c_h0_count)

			if haplo_count_ratio > max_haplo_count_ratio:
				# the haplotype of ss reads is not distinguishable enough based in the ss map dir
				continue

			w_haplo = 0 if w_h0_c_h1_count < w_h1_c_h0_count else 1
			c_haplo = 1 - w_haplo

			# TODO: replace it with assert statement to make sure every chrom is clustered in both directions
			if chrom + '_16' not in chrom_dir_to_clust or chrom + '_0' not in chrom_dir_to_clust:
				print('Warning: chrom ' + chrom + ' is not properly clustered...')
				continue

			chrom_w_clust = chrom_dir_to_clust[chrom + '_16']
			chrom_c_clust = chrom_dir_to_clust[chrom + '_0']

			print(chrom_w_clust + '\t' + str(w_haplo) + '\t' + chrom + '_16', file=out)
			print(chrom_c_clust + '\t' + str(c_haplo) + '\t' + chrom + '_0' , file=out)


