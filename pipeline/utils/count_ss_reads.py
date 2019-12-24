import pysam

def output_ss_counts_in_chroms(ss_haplotagged_bam, ss_counts_file):

	valid_chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX']
	ss_counts_in_chroms = {chrom:{'w':0, 'c':0, 'w_h0':0, 'w_h1':0, 'c_h0':0, 'c_h1':0} for chrom in valid_chroms}

	bamfile = pysam.AlignmentFile(ss_haplotagged_bam, 'rb')

	for read in bamfile.fetch():
		chrom = read.reference_name

		if chrom not in valid_chroms:
			continue

		read_class = []

		dir = 'w' if read.is_reverse else 'c'

		ss_counts_in_chroms[chrom][dir] += 1

		if read.has_tag("HP"):
			hp = str(read.get_tag("HP")-1)
			ss_counts_in_chroms[chrom][dir + '_h' + hp] += 1

	
	with open(ss_counts_file, 'w') as out:
		print('chrom\tW\tC\tW_h1\tW_h2\tC_h1\tC_h2', file=out)

		for chrom in ss_counts_in_chroms:
			ss_chrom_counts = ss_counts_in_chroms[chrom]
			counts = [str(ss_chrom_counts[ss_class]) for ss_class in ['w','c','w_h0','w_h1','c_h0','c_h1']]

			print('\t'.join([chrom] + counts), file=out)

	return ss_counts_in_chroms



