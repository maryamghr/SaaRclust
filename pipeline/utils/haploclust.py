from __future__ import division
from bubble_long_read_alignment import *
from string_functions import *
from parsing import *
import pysam
import numpy
import time
import pdb
from argparse import ArgumentParser
import gzip

#bubble_10_allele_1_chain_1_km_10.5_node_13237557        125     0       125     +       m54329U_190607_185248/15794285/ccs_0_chr4_10075151      10886   10605   10730   124
#     125     21      NM:i:1  ms:i:244        AS:i:244        nn:i:0  tp:A:P  cm:i:81 s1:i:124        s2:i:108        de:f:0.0080     cg:Z:125M
#1	string	Query sequence name
#2	int	Query sequence length
#3	int	Query start coordinate (0-based)
#4	int	Query end coordinate (0-based)
#5	char	‘+’ if query/target on the same strand; ‘-’ if opposite
#6	string	Target sequence name
#7	int	Target sequence length
#8	int	Target start coordinate on the original strand
#9	int	Target end coordinate on the original strand
#10	int	Number of matching bases in the mapping
#11	int	Number bases, including gaps, in the mapping
#12	int	Mapping quality (0-255 with 255 for missing)


def add_bubbles_het_positions(bubbles):
	for bubble_id, bubble in bubbles.items():
		bubble.add_het_positions()
		
#def reset_num_haplo_long_reads(bubbles):
#	for bubble_id, bubble in bubbles.items():
#		bubble.allele0.num_haplo_long_reads = [0,0]
#		bubble.allele1.num_haplo_long_reads = [0,0]

def haploclust_long_reads(long_reads):
	for read_name, long_read in long_reads.items():
		long_read.phase()
		
def haploclust_bubbles(bubbles, test_bubbles):
	for bubble_id, bubble in bubbles.items():
		#if bubble_id in test_bubbles:
		#	pdb.set_trace()
		bubble.phase()
	
def output_kmers(long_reads, kmers_file):
	print('outputting the kmers')
	with open(kmers_file, 'w') as out:
		for read_name, long_read in long_reads.items():
			for bubble_allele, aln in long_read.alignments.items():
				print(aln.output_kmers(), file=out)
	
def iterative_haplo_clust(bubbles, long_reads, q, itr=2):
	
	'''
	Assign haplotypes to long reads and bubbles in {itr} iterations, starting from the first set of phased bubbles from the input file
	
	Parameters:
		itr (int)						: the number of halotype clustering iterations
		bubble_first_itr_phase_file: Path to the file containing the first set of phased bubbles
		bubbles			 				: A dictionary {bubble_id -> bubbles}
		long_reads		 				: A dictionary {long_read_name -> long_read}
	'''
	test_bubbles = {}
	
	#add_bubble_allele_pred_haplo(bubble_first_itr_phase_file, bubbles)
	add_bubbles_het_positions(bubbles)
	
	print('haploclust long reads first iteration...')
	for read_name, long_read in long_reads.items():
		long_read.set_alignments_edit_dist(q)
		long_read.phase()
		
		# testing:
		if long_read.pred_haplo != None:
			for bubble_allele, aln in long_read.alignments.items():
				if bubble_allele.pred_haplo == None and bubble_allele.id == 0:
					test_bubbles[bubble_allele.bubble.id] = bubble_allele.bubble
	
	for it in range(itr-1):
		print('haploclust iteration', it+2)
		haploclust_bubbles(bubbles, test_bubbles)
		#reset_num_haplo_long_reads(bubbles)
		haploclust_long_reads(long_reads)
		
	return test_bubbles

def output_phasing(bubbles, long_reads, test_bubbles, bubble_phasing_file, long_read_phasing_file, test_bubble_phasing_file):
	with open(bubble_phasing_file, 'w') as out:
		print("bubble_id\tallele0haplo", file=out)
		for bubble_id, bubble in bubbles.items():
			print(bubble_id, "\t", bubble.allele0.pred_haplo, file=out)

	with open(long_read_phasing_file, 'w') as out:
		print("long_read_name\thaplo0_edit_dist\thaplo1_edit_dist\thaplotype", file=out)
		for read_name, long_read in long_reads.items():
			print(read_name, "\t", long_read.haplo0_edit_dist, "\t", long_read.haplo1_edit_dist, "\t", long_read.pred_haplo, file=out)
			
	with open(test_bubble_phasing_file, 'w') as out:
		print("bubble_id\tallele0haplo", file=out)
		for bubble_id, bubble in test_bubbles.items():
			print(bubble_id, "\t", bubble.allele0.pred_haplo, file=out)
