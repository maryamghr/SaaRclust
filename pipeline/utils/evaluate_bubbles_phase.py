from __future__ import division
import pysam

global valid_chroms
valid_chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX']

def print_dict_head(dic, num):
	for i in range(num):
		key = list(dic.keys())[i]
		print(key, ':', dic[key])

def get_ground_true_bubbles_chrom_and_h0allele(bubble_haplotagged_bam_file):
	num_bubbles = 0
	bubble_id_to_chrom, bubble_id_to_h0_allele = {}, {}
	samfile = pysam.AlignmentFile(bubble_haplotagged_bam_file, 'rb')
	
	for read in samfile.fetch():
		
		bubble_name = read.query_name
		bubble_name_sp = bubble_name.split('_')
		bubble_id, bubble_allele = int(bubble_name_sp[1]), int(bubble_name_sp[3])-1

		if bubble_allele==0:
			num_bubbles += 1

		if not read.has_tag("HP"):
			# the read is not haplotagged
			continue

		chrom, haplotype = read.reference_name, read.get_tag("HP")-1
		
		if not chrom.startswith('chr'):
			# the chromosome name is not valid
			continue
		
		bubble_h0_allele = bubble_allele if haplotype==0 else 1-bubble_allele

		if bubble_id not in bubble_id_to_chrom:
			bubble_id_to_chrom[bubble_id] = chrom
			bubble_id_to_h0_allele[bubble_id] = bubble_h0_allele

	return num_bubbles, bubble_id_to_chrom, bubble_id_to_h0_allele


def get_clust_to_chrom(clust_to_chrom_file):
	clust_to_chrom = {}
	with open(clust_to_chrom_file) as f:
		for line in f:
			if not line.startswith('chr'):
				# it is a hedear line
				continue

			chrom, clust_forward, clust_backward = line.split()

			if clust_forward not in clust_to_chrom:
				# remove the direction from chrom name
				clust_to_chrom[clust_forward] = chrom.split('_')[0]

	return clust_to_chrom

def get_bubble_id_to_clust(bubble_clust_file):
	bubble_id_to_clust = {}
	with open(bubble_clust_file) as f:
		for line in f:
			sp = line.split()

			if not sp[0].startswith('V'):
				# it is a header line
				continue

			bubble_id_to_clust[int(sp[1])] = sp[0]

	return bubble_id_to_clust


def get_bubble_id_to_h0_allele(bubble_phase_file):
	bubble_id_to_h0_allele = {}
	with open(bubble_phase_file) as f:
		for line in f:
			# skip the header line
			if line.startswith('b'):
				continue

			sp = line.split()

			bubble_id_to_h0_allele[int(sp[0])] = int(sp[1])

	return bubble_id_to_h0_allele



def evaluate_bubble_phase(num_bubbles, bubble_id_to_chrom, ground_true_bubble_id_to_h0_allele, clust_to_chrom, bubble_id_to_clust, bubble_id_to_h0_allele, output_file, false_phased_bubbles_file):
	
	num_clustered_bubbles = 0 #len(bubble_id_to_clust)
	num_true_clustered_bubbles, num_phased_bubbles, num_true_phased_bubbles = 0, 0, 0
	chrom_to_num_bubbles_with_same_haplo = {}
	chrom_to_num_bubbles_with_diff_haplo = {}
	chrom_to_flip_haplotype = {}
	chrom_to_num_clustered_bubbles = {chrom:0 for chrom in valid_chroms}
	chrom_to_num_true_clustered_bubbles = {chrom:0 for chrom in valid_chroms}
	chrom_to_num_phased_bubbles = {chrom:0 for chrom in valid_chroms}
	chrom_to_num_true_phased_bubbles = {chrom:0 for chrom in valid_chroms}

	

	for bubble_id in bubble_id_to_clust:
		clust = bubble_id_to_clust[bubble_id]


		if clust not in clust_to_chrom or bubble_id not in bubble_id_to_chrom:
			continue

		ground_true_chrom = bubble_id_to_chrom[bubble_id]

		num_clustered_bubbles += 1
		chrom_to_num_clustered_bubbles[ground_true_chrom] += 1

		if clust_to_chrom[clust] != ground_true_chrom:
			continue

		num_true_clustered_bubbles += 1
		chrom_to_num_true_clustered_bubbles[ground_true_chrom] += 1

		if bubble_id not in bubble_id_to_h0_allele:
			# bubble is not phased
			continue

		num_phased_bubbles += 1
		chrom_to_num_phased_bubbles[ground_true_chrom] += 1

		if ground_true_chrom not in chrom_to_num_bubbles_with_same_haplo:
			chrom_to_num_bubbles_with_same_haplo[ground_true_chrom], chrom_to_num_bubbles_with_diff_haplo[ground_true_chrom] = 0, 0

		if bubble_id not in ground_true_bubble_id_to_h0_allele:
			continue
			
		if bubble_id_to_h0_allele[bubble_id] == ground_true_bubble_id_to_h0_allele[bubble_id]:
			chrom_to_num_bubbles_with_same_haplo[ground_true_chrom] += 1
		else:
			chrom_to_num_bubbles_with_diff_haplo[ground_true_chrom] += 1



	for chrom in chrom_to_num_bubbles_with_same_haplo:
		chrom_to_num_true_phased_bubbles[chrom] = max(chrom_to_num_bubbles_with_same_haplo[chrom], chrom_to_num_bubbles_with_diff_haplo[chrom])
		num_true_phased_bubbles += chrom_to_num_true_phased_bubbles[chrom]

		chrom_to_flip_haplotype[chrom] = False if chrom_to_num_bubbles_with_same_haplo[chrom] > chrom_to_num_bubbles_with_diff_haplo[chrom] else True

	with open(false_phased_bubbles_file, 'w') as out_false:
		
		print("bubble_id\th0_allele\tground_true_h0_allele_after_possibly_flipping_haplo\tground_true_id_to_h0_allele\tchrom_to_flip_haplotype\tchrom", file=out_false)
		
		for bubble_id in bubble_id_to_h0_allele:
			
			if bubble_id not in bubble_id_to_clust:
				continue

			clust = bubble_id_to_clust[bubble_id]

			if clust not in clust_to_chrom:
				continue

			if bubble_id not in ground_true_bubble_id_to_h0_allele:
				continue

			chrom = clust_to_chrom[clust]

			ground_true_h0_allele = ground_true_bubble_id_to_h0_allele[bubble_id] if not chrom_to_flip_haplotype[chrom] else 1-ground_true_bubble_id_to_h0_allele[bubble_id]

			if bubble_id_to_h0_allele[bubble_id] != ground_true_h0_allele:
				print(str(bubble_id) + "\t" + str(bubble_id_to_h0_allele[bubble_id]) + "\t" + str(ground_true_h0_allele) + "\t" + str(ground_true_bubble_id_to_h0_allele[bubble_id]) + "\t" + str(chrom_to_flip_haplotype[chrom]) + "\t" + chrom, file=out_false)

		

	print('chrom_to_num_bubbles_with_same_haplo')
	print(chrom_to_num_bubbles_with_same_haplo)
	print('chrom_to_num_bubbles_with_diff_haplo')
	print(chrom_to_num_bubbles_with_diff_haplo)		

	with open(output_file, 'w') as out:
		print('total number of bubbles = ' + str(num_bubbles), file=out)
		print('number of clustered bubbles = ' + str(num_clustered_bubbles) + ', (' + str(num_clustered_bubbles*100/num_bubbles) + ' % of #bubbles)', file=out)
		print('number of true clustered bubbles = ' + str(num_true_clustered_bubbles) + ', (' + str(num_true_clustered_bubbles*100/num_clustered_bubbles) + ' % of #clustered bubbled)', file=out)
		print('number of phased bubbles = ' + str(num_phased_bubbles) + ', (' + str(num_phased_bubbles*100/num_bubbles) + '% of #bubbles)', file=out)
		print('number of true phased bubbles = ' + str(num_true_phased_bubbles) + ', (' + str(num_true_phased_bubbles*100/num_phased_bubbles) + '% of #phased bubbles)', file=out)
		print('chromosome wise clustering accuracy:')
		print('chrom\t#clustered_bubbles\ttrue_clustered_bubbles\tclustering_accuracy\tnum_phased_bubbles\tnum_true_phased_bubbles\tphasing_accuracy', file=out)
		for chrom in chrom_to_num_true_clustered_bubbles:
			num_clustered, num_true_clustered, num_phased, num_true_phased = chrom_to_num_clustered_bubbles[chrom], chrom_to_num_true_clustered_bubbles[chrom], chrom_to_num_phased_bubbles[chrom], chrom_to_num_true_phased_bubbles[chrom]
			print(chrom, '\t', num_clustered, '\t', num_true_clustered, '\t', num_true_clustered*100/num_clustered, num_phased, '\t', num_true_phased, '\t', num_true_phased*100/num_phased, file=out)




