from __future__ import division
import pysam


###############################################
#### evaluating long reads phasing info #######
###############################################

global valid_chroms
valid_chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX']

def get_ground_true_chrom_haplo(long_reads_haplotagged_bam_files):
	chrom_to_num_ground_true_phased_reads = {chrom:0 for chrom in valid_chroms}
	ground_true_read_to_chrom, ground_true_read_to_haplo = {}, {}
	
	for alignmentfile in long_reads_haplotagged_bam_files:
		print('processing', alignmentfile, '...')
		samfile = pysam.AlignmentFile(alignmentfile, 'rb')
		for read in samfile.fetch():
		
			read_name = read.query_name
		
			if not read.has_tag("HP"):
				# the read is not haplotagged
				continue

			chrom, haplotype = read.reference_name, read.get_tag("HP")-1
		
			if not chrom in valid_chroms:
				# the chromosome name is not valid
				continue
		
			chrom_to_num_ground_true_phased_reads[chrom] += 1
			ground_true_read_to_chrom[read_name], ground_true_read_to_haplo[read_name] = chrom, haplotype

	return chrom_to_num_ground_true_phased_reads, ground_true_read_to_chrom, ground_true_read_to_haplo


def get_reads_haplotypes(long_reads_phase_file):
	read_to_haplo = {}
	with open(long_reads_phase_file) as f:
		# skip the header line
		next(f)
		for line in f:
			sp = line.split()
			read_to_haplo[sp[0]] = sp[-1]

	return read_to_haplo


def output_long_reads_phasing_accuracy(chrom_to_num_ground_true_phased_reads, ground_true_read_to_chrom, ground_true_read_to_haplo,  read_to_haplo, output_file):
	
	chrom_to_num_predicted_phased_reads = {chrom:0 for chrom in valid_chroms}
	chrom_to_num_correctly_phased_reads = {chrom:0 for chrom in valid_chroms}

	for read in read_to_haplo:
		if read not in ground_true_read_to_haplo:
			# the read is not phased in the ground true data
			continue

		chrom = ground_true_read_to_chrom[read]
		actual_haplo = ground_true_read_to_haplo[read]
		predicted_haplo = read_to_haplo[read]

		if read_to_haplo[read] == "?":
			# the read is not phased
			continue

		chrom_to_num_predicted_phased_reads[chrom] += 1

		if int(predicted_haplo) == actual_haplo:
			chrom_to_num_correctly_phased_reads[chrom] += 1

	
	with open(output_file, 'w') as out:
		print("chrom\thaplo_score\tswitch_haplo_score\tnum_actual_phased_reads\tnum_predicted_phased_reads\tnum_correctly_phased_reads\tfraction_phased_reads\tfraction_correctly_phased_reads", file=out)
		for chrom in valid_chroms:
			# change the score (#correctly phased reads) if switching h0 and h1 leads to a better score (the order of h0 and h1 does not matter)
			haplo_score, switch_haplo_score = chrom_to_num_correctly_phased_reads[chrom], chrom_to_num_predicted_phased_reads[chrom]-chrom_to_num_correctly_phased_reads[chrom]
			num_correctly_phase_reads = max(haplo_score, switch_haplo_score)
			num_actual_phased_reads = chrom_to_num_ground_true_phased_reads[chrom]
			num_predicted_phased_reads = chrom_to_num_predicted_phased_reads[chrom]
			fraction_phased_reads = chrom_to_num_predicted_phased_reads[chrom]/chrom_to_num_ground_true_phased_reads[chrom]
			fraction_correctly_phased_reads = num_correctly_phase_reads / chrom_to_num_predicted_phased_reads[chrom]

			print(chrom, '\t', haplo_score, '\t', switch_haplo_score, '\t', num_actual_phased_reads, '\t', num_predicted_phased_reads, '\t', num_correctly_phase_reads, '\t', fraction_phased_reads, '\t', fraction_correctly_phased_reads, file=out)
			
			
			


