from __future__ import division
from bubble_long_read_alignment import *
from parsing import *
#import pysam
#import numpy
#import time
#import pdb
from argparse import ArgumentParser

global valid_chroms
valid_chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX']


def evaluate_bubble_clustering(bubbles, clust_to_chrom, output_file):

	'''
	'''

	start_time = time.time()
	print('evaluating bubbles clustering')

	num_bubbles=len(bubbles)
	num_chr_clustered_bubbles=0
	num_garbage_chr_clustered_bubbles=0
	num_true_chr_clustered_bubbles=0
	num_haplo_clustered_bubbles=0
	num_true_haplo_clustered_bubbles=0
	num_haploclust_false_pos=0
	num_haploclust_false_neg=0

	# defining the same values per chromosome
	chrom_num_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_chr_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_garbage_chr_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_true_chr_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_haplo_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_true_haplo_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_false_haplo_clustered_bubbles={chrom:0 for chrom in valid_chroms}
	chrom_num_haploclust_false_pos={chrom:0 for chrom in valid_chroms}
	chrom_num_haploclust_false_neg={chrom:0 for chrom in valid_chroms}
	chrom_switch_haplo={chrom:False for chrom in valid_chroms}
	
	#for bubble_id, bubble in bubbles.items():
	for bubble_id in bubbles:
		bubble = bubbles[bubble_id]

		if bubble.actual_chrom == None or bubble.actual_chrom not in valid_chroms:
			# the bubble does not have a valid actual chrom
			continue

		chrom = bubble.actual_chrom

		chrom_num_bubbles[chrom]+=1

		if bubble.clust == None:
			bubble.pred_type="not_chrom_clust"
			continue

		if bubble.clust not in clust_to_chrom:
			num_garbage_chr_clustered_bubbles+=1
			chrom_num_garbage_chr_clustered_bubbles[chrom]+=1
			bubble.pred_type="garbage_clust"
			continue

		num_chr_clustered_bubbles+=1
		chrom_num_chr_clustered_bubbles[chrom]+=1

		if clust_to_chrom[bubble.clust] == bubble.actual_chrom:
			num_true_chr_clustered_bubbles+=1
			chrom_num_true_chr_clustered_bubbles[chrom]+=1

		if bubble.allele0.pred_haplo == None:
			bubble.pred_type="not_haplo_clust"

			if bubble.allele0.actual_haplo != None:
				bubble.pred_type="haplo_clust_false_neg"
				num_haploclust_false_neg+=1
				chrom_num_haploclust_false_neg[chrom]+=1

			continue

		num_haplo_clustered_bubbles+=1
		chrom_num_haplo_clustered_bubbles[chrom]+=1

		if bubble.allele0.actual_haplo==None:
			bubble.pred_type="haploclust_false_pos"
			num_haploclust_false_pos+=1
			chrom_num_haploclust_false_pos[chrom]+=1
			continue

		elif bubble.allele0.pred_haplo == bubble.allele0.actual_haplo:
			bubble.pred_type="true_haplo_clust"
			chrom_num_true_haplo_clustered_bubbles[chrom]+=1
			continue

		else:
			bubble.pred_type="false_haplo_clust"
			chrom_num_false_haplo_clustered_bubbles[chrom]+=1

	# revise chrom num true haplotypes after switching the haplotypes in chroms if necessary
	for chrom in valid_chroms:
		true_haplo_clust=chrom_num_true_haplo_clustered_bubbles[chrom]
		false_haplo_clust=chrom_num_false_haplo_clustered_bubbles[chrom]
		
		if true_haplo_clust < false_haplo_clust:
			# switch haplotypes in the chrom
			chrom_switch_haplo[chrom]=True
			true_haplo_clust, false_haplo_clust = false_haplo_clust, true_haplo_clust
			chrom_num_true_haplo_clustered_bubbles[chrom]=true_haplo_clust
			chrom_num_false_haplo_clustered_bubbles[chrom]=false_haplo_clust
			
		num_true_haplo_clustered_bubbles+=true_haplo_clust

	# revise the type of the bubble if the haplotype is switched in the chromosome
	for bubble_id, bubble in bubbles.items():
		if bubble.pred_type=="true_haplo_clust" or bubble.pred_type=="false_haplo_clust":
			if chrom_switch_haplo[bubble.actual_chrom]:
				bubble.pred_type="true_haplo_clust" if bubble.pred_type=="false_haplo_clust" else "false_haplo_clust"

	# writing the performance statistics in the outout file
	with open(output_file, 'w') as out:
		print('*** Note: false positive haplotype clustered bubbles are also counted in haplotype clustering accuracy', file=out)
		print('total number of bubbles =', num_bubbles, file=out)
		print('number of chrom clustered bubbles =', num_chr_clustered_bubbles, ', (', num_chr_clustered_bubbles*100/num_bubbles, ' % of #bubbles)', file=out)
		print('chrom clustering accuracy =', num_true_chr_clustered_bubbles*100/num_chr_clustered_bubbles, file=out)
		print('number of haplo clustered bubbles = ', num_haplo_clustered_bubbles, ', (', num_haplo_clustered_bubbles*100/num_bubbles, '% of #bubbles)', file=out)
		print('number of true haplo clustered bubbles = ', num_true_haplo_clustered_bubbles, file=out)
		print('haplo clustering accuracy =', num_true_haplo_clustered_bubbles*100/num_haplo_clustered_bubbles, file=out)
		print('haplo clustering false positive rate =', num_haploclust_false_pos*100/num_bubbles, file=out)
		print('haplo clustering false negative rate =', num_haploclust_false_neg*100/num_bubbles, file=out)


		print('chromosome wise clustering accuracy:', file=out)

		print('chrom\t#bubbles\t\
		#chr_clustered_bubbles\tfraction_chr_clustered_bubbles\tchr_clustering_accuracy\t\
		#haplo_clustered_bubbles\tfraction_haplo_clustered_bubbles\thaplo_clustering_accuracy\t\
		false_pos\tfalse_neg', file=out)

		for chrom in valid_chroms:
			print(chrom, '\t', chrom_num_bubbles[chrom], '\t', 
			chrom_num_chr_clustered_bubbles[chrom], '\t', 
			chrom_num_chr_clustered_bubbles[chrom]*100/chrom_num_bubbles[chrom], '\t', 
			chrom_num_true_chr_clustered_bubbles[chrom]*100/chrom_num_chr_clustered_bubbles[chrom], '\t',
			chrom_num_haplo_clustered_bubbles[chrom], '\t', 
			chrom_num_haplo_clustered_bubbles[chrom]*100/chrom_num_bubbles[chrom], '\t', 
			chrom_num_true_haplo_clustered_bubbles[chrom]*100/chrom_num_haplo_clustered_bubbles[chrom], '\t', 
			chrom_num_haploclust_false_pos[chrom]*100/chrom_num_bubbles[chrom], '\t', 
			chrom_num_haploclust_false_neg[chrom]*100/chrom_num_bubbles[chrom], '\t', 		
			file=out)

	print('elapsed time =', time.time()-start_time)	


def output_bubbles_haplo_dist(bubbles, output_file, with_km=True):

# TODO: assert that alignments are present
	'''
	'''

	start_time = time.time()
	print('outputting the bubbles haplo edit dist')

	with open(output_file, 'w') as out:
		km = "km" if with_km else ""
		print('chrom\tbubble_id\tdist_h0\tdist_h1\tnum_al0_h0_long_reads\tnum_al0_h1_long_reads\tactual_type\tpred_type\t'+str(km), file=out)

		for bubble_id, bubble in bubbles.items():

			h0_edit_dist, h1_edit_dist = bubble.allele0.get_haplotypes_edit_dist()		
			#h0_num_aln_reads, h1_num_aln_reads = bubble.get_haplo_num_aligned_reads()

			km = bubble.allele0.km+bubble.allele1.km if with_km else ""

			print(str(bubble.actual_chrom), '\t', bubble_id, '\t', 
						h0_edit_dist, '\t', h1_edit_dist, '\t', 
						bubble.num_al0_h0_reads, '\t', bubble.num_al0_h1_reads, '\t', 
						bubble.actual_type, '\t', bubble.pred_type, '\t', km, file=out)

	print('elapsed time =', time.time()-start_time)


def evaluate_long_read_clustering(long_reads, output_file):

	'''
	'''

	start_time = time.time()
	print('evaluating long reads clustering')

	num_long_reads=len(long_reads)
	num_haplo_clustered_long_reads=0
	num_true_haplo_clustered_long_reads=0
	num_haploclust_false_pos=0
	num_haploclust_false_neg=0

	# defining the same values per chromosome
	chrom_num_long_reads={chrom:0 for chrom in valid_chroms}
	chrom_num_haplo_clustered_long_reads={chrom:0 for chrom in valid_chroms}
	chrom_num_true_haplo_clustered_long_reads={chrom:0 for chrom in valid_chroms}
	chrom_num_false_haplo_clustered_long_reads={chrom:0 for chrom in valid_chroms}
	chrom_num_haploclust_false_pos={chrom:0 for chrom in valid_chroms}
	chrom_num_haploclust_false_neg={chrom:0 for chrom in valid_chroms}
	chrom_switch_haplo={chrom:False for chrom in valid_chroms}
	
	for read_name, long_read in long_reads.items():

		if long_read.actual_chrom == None or long_read.actual_chrom not in valid_chroms:
			# the long_read does not have a valid actual chrom
			continue

		chrom = long_read.actual_chrom

		chrom_num_long_reads[chrom]+=1

		if long_read.pred_haplo == None:
			long_read.pred_type="not_haplo_clust"

			if long_read.actual_haplo != None:
				long_read.pred_type="haplo_clust_false_neg"
				num_haploclust_false_neg+=1
				chrom_num_haploclust_false_neg[chrom]+=1

			continue

		num_haplo_clustered_long_reads+=1
		chrom_num_haplo_clustered_long_reads[chrom]+=1

		if long_read.actual_haplo==None:
			long_read.pred_type="haploclust_false_pos"
			num_haploclust_false_pos+=1
			chrom_num_haploclust_false_pos[chrom]+=1
			continue

		elif long_read.pred_haplo == long_read.actual_haplo:
			long_read.pred_type="true_haplo_clust"
			chrom_num_true_haplo_clustered_long_reads[chrom]+=1
			continue

		else:
			long_read.pred_type="false_haplo_clust"
			chrom_num_false_haplo_clustered_long_reads[chrom]+=1

	# revise chrom num true haplotypes after switching the haplotypes in chroms if necessary
	for chrom in valid_chroms:
		true_haplo_clust=chrom_num_true_haplo_clustered_long_reads[chrom]
		false_haplo_clust=chrom_num_false_haplo_clustered_long_reads[chrom]
		
		if true_haplo_clust < false_haplo_clust:
			# switch haplotypes in the chrom
			chrom_switch_haplo[chrom]=True
			true_haplo_clust, false_haplo_clust = false_haplo_clust, true_haplo_clust
			chrom_num_true_haplo_clustered_long_reads[chrom]=true_haplo_clust
			chrom_num_false_haplo_clustered_long_reads[chrom]=false_haplo_clust
			
		num_true_haplo_clustered_long_reads+=true_haplo_clust
		
	print('switch haplotypes:')
	print(chrom_switch_haplo)

	# revise the type of the long_read if the haplotype is switched in the chromosome
	for read_name, long_read in long_reads.items():
		if long_read.pred_type=="true_haplo_clust" or long_read.pred_type=="false_haplo_clust":
			if chrom_switch_haplo[long_read.actual_chrom]:
				long_read.pred_type="true_haplo_clust" if long_read.pred_type=="false_haplo_clust" else "false_haplo_clust"

	# writing the performance statistics in the outout file
	with open(output_file, 'w') as out:
		print('*** Note: false positive haplotype clustered long_reads are also counted in haplotype clustering accuracy', file=out)
		print('total number of long_reads =', num_long_reads, file=out)
		print('number of haplo clustered long_reads = ', num_haplo_clustered_long_reads, ', (', num_haplo_clustered_long_reads*100/num_long_reads, '% of #long_reads)', file=out)
		print('number of true haplo clustered =', num_true_haplo_clustered_long_reads, file=out)
		print('haplo clustering accuracy =', num_true_haplo_clustered_long_reads*100/num_haplo_clustered_long_reads, file=out)
		print('haplo clustering false positive rate =', num_haploclust_false_pos*100/num_long_reads, file=out)
		print('haplo clustering false negative rate =', num_haploclust_false_neg*100/num_long_reads, file=out)


		print('chromosome wise clustering accuracy:')

		print('chrom\t#long_reads\t\
		#haplo_clustered_long_reads\tfraction_haplo_clustered_long_reads\thaplo_clustering_accuracy\t\
		false_pos\tfalse_neg', file=out)

		for chrom in valid_chroms:
			print(chrom, '\t', chrom_num_long_reads[chrom], '\t', 
			#chrom_num_chr_clustered_long_reads[chrom], '\t', 
			#chrom_num_chr_clustered_long_reads[chrom]*100/chrom_num_long_reads[chrom], '\t', 
			#chrom_num_true_chr_clustered_long_reads[chrom]*100/chrom_num_chr_clustered_long_reads[chrom], '\t',
			chrom_num_haplo_clustered_long_reads[chrom], '\t', 
			chrom_num_haplo_clustered_long_reads[chrom]*100/chrom_num_long_reads[chrom], '\t', 
			chrom_num_true_haplo_clustered_long_reads[chrom]*100/chrom_num_haplo_clustered_long_reads[chrom], '\t', 
			chrom_num_haploclust_false_pos[chrom]*100/chrom_num_long_reads[chrom], '\t', 
			chrom_num_haploclust_false_neg[chrom]*100/chrom_num_long_reads[chrom], '\t', 		
			file=out)

	print('elapsed time =', time.time()-start_time)	


def output_long_reads_haplo_dist(long_reads, output_file):

	# TODO: assert that alignments are present
	'''
	'''

	start_time = time.time()
	print('outputting the long reads haplo edit dist')

	with open(output_file, 'w') as out:
		print('chrom\tread_name\tdist_h0\tdist_h1\tnum_h0_bubbles\tnum_h1_bubbles\tactual_type\tpred_type\tpred_haplo', file=out)

		for read_name, long_read in long_reads.items():

			long_read.set_haplotypes_edit_dist()
			num_aln_reads = len(long_read.alignments)

			print(str(long_read.actual_chrom), '\t', read_name, '\t', 
						long_read.haplo0_edit_dist, '\t', long_read.haplo1_edit_dist, '\t', 
						long_read.num_haplo_bubbles[0], '\t', long_read.num_haplo_bubbles[1], '\t',
						long_read.actual_type, '\t', long_read.pred_type, long_read.pred_haplo, file=out)

	print('elapsed time =', time.time()-start_time)
	

def output_sampled_long_reads(num_sample, edit_dist_fraction_range, file_name):
	n = 0
	
	with open(file_name, 'w') as out:
		print("bubbleName\tbubbleAllele\tPBname\tbubbleKmer\tPBkmer\tkmersEditDistance\tbubble_alle_pred_haplo\tlong_read_pred_haplo", file=out)
		
		for read_name, long_read in long_reads.items():
			if n > num_sample:
					break
					
			d0, d1 = long_read.haplo0_edit_dist, long_read.haplo1_edit_dist
			if d0+d1 == 0:
				continue
			
			dist_frac = min(d0, d1)/(d0+d1)
		
			if edit_dist_fraction_range[0] < dist_frac < edit_dist_fraction_range[1]:
				n += 1
				
				for aln in long_read.alignments:
						print(aln.output_print(), file=out)
						
				print('h_dist0 =', d0, file=out)
				print('h_dist1 =', d1, '\n', file=out)
						

if __name__ == "__main__":
	parser = ArgumentParser(description=__doc__)
	parser.add_argument("--bubble_bam_file", type=str, help="Bubble haplotagged bam file", required=True)
	parser.add_argument("--bubble_phase_file", type=str, help="Bubbles haplo clustering file", required=True)
	parser.add_argument("--bubble_clust_file", type=str, help="Bubbles cluster file", required=True)
	parser.add_argument("--clust_to_chrom_file", type=str, help="Cluster to chrom mapping file", required=True)
	parser.add_argument("--het_kmers_files", nargs='*', help="The set of bubble/long reads het kmers files", required=True)
	parser.add_argument("--long_read_haplotagged_bam_files", nargs='*', help="The set of long reads haplotagged bam files", required=True)
	parser.add_argument("--long_read_phase_files", nargs='*', help="The set of long reads haplo clustering files", required=True)
	parser.add_argument("--bubbles_clusering_evaluation_file", type=str, help="The output bubbles clustring evaluation file", required=True)
	parser.add_argument("--bubbles_haplo_edit_dist_file", type=str, help="The output bubbles haplo edit dist file", required=True)
	parser.add_argument("--long_read_phase_evaluation_file", type=str, help="The output long reads phasing evaluation file", required=True)
	parser.add_argument("--long_reads_haplo_edit_dist_file", type=str, help="The output long reads haplo edit dist file", required=True)

	# testing output files for observation 
	parser.add_argument("--long_reads_with_small_frac_haplo_edit_dist", type=str, help="The output sample long reads with small fraction of haplo_edit_dist", required=True)
	parser.add_argument("--long_reads_with_peak_frac_haplo_edit_dist", type=str, help="The output sample long reads with peak fraction of haplo_edit_dist", required=True)
	######################################
	
	parser.add_argument("--with_km", action='store_true', help="True if bubbles have km information in their name")
	args = parser.parse_args()

	#counts(args.sample, args.input_bam, args.input_bed, args.counts_output)
	with_km = True if "with_km" in args else False
	print('with_km', 	with_km)

	bubbles = get_bubbles_from_bam(args.bubble_bam_file, with_km)
	clust_to_chrom = get_clust_to_chrom(args.clust_to_chrom_file)
	add_bubble_clust(args.bubble_clust_file, bubbles)
	add_bubble_allele_pred_haplo(args.bubble_phase_file, bubbles)
	long_reads = get_long_reads_from_bam(args.long_read_haplotagged_bam_files)
	add_long_reads_pred_haplotype(args.long_read_phase_files, long_reads)
	set_alignments_from_kmers_file(args.het_kmers_files, bubbles, long_reads)
	evaluate_bubble_clustering(bubbles, clust_to_chrom, args.bubbles_clusering_evaluation_file)
	output_bubbles_haplo_dist(bubbles, args.bubbles_haplo_edit_dist_file, with_km)
	evaluate_long_read_clustering(long_reads, args.long_read_phase_evaluation_file)
	output_long_reads_haplo_dist(long_reads, args.long_reads_haplo_edit_dist_file)
	
	# sampling long reads:
	#output_sampled_long_reads(100, (0.35, 0.4), args.long_reads_with_peak_frac_haplo_edit_dist)
	#output_sampled_long_reads(10, (0, 0.1), args.long_reads_with_small_frac_haplo_edit_dist)



