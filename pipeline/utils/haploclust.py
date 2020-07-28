from __future__ import division
from bubble_long_read_alignment import *
from string_functions import *
from parsing import *
from evaluate_haploclust import *
import pysam
import numpy
import time
import pdb
from argparse import ArgumentParser
import gzip
import math

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
		

def haploclust_long_reads(long_reads, min_haplotagged_bubbles=1):
	start_time = time.time()
	num = 0
	total_reads = len(long_reads)-(len(long_reads)%10)
	for read_name, long_read in long_reads.items():
		num += 1
		if num*10 % total_reads == 0:
			print(num*100/total_reads, '% of reads processed')
		long_read.phase(min_haplotagged_bubbles=min_haplotagged_bubbles)
		
	print('elapsed time =', time.time()-start_time, 's')

	
def haploclust_bubbles(bubbles):
	start_time = time.time()
	num = 0
	total_bubbles = len(bubbles)-(len(bubbles)%10)
	for bubble_id, bubble in bubbles.items():
		num += 1
		if num*10 % total_bubbles == 0:
			print(num*100/total_bubbles, '% of bubbles processed')
		bubble.phase()
		
	print('elapsed time =', time.time()-start_time, 's')


def output_kmers(long_reads, kmers_file):
	print('outputting the kmers')
	num = 0
	total_reads = len(long_reads)-(len(long_reads)%10)
	with open(kmers_file, 'w') as out:
		for read_name, long_read in long_reads.items():
			num += 1
			if num*10 % total_reads == 0:
				print(num*100/total_reads, '% of reads processed')
			
			for aln in long_read.alignments:
				print(aln.output_kmers(), file=out)

				
def get_km_quantile(bubbles,p):
	km = []
	for bubble_id, bubble in bubbles.items():
		km.append(bubble.allele0.km)
		km.append(bubble.allele1.km)
		
	km.sort()
	
	idx = math.floor(p*len(km))
	
	return km[idx]
	

def trim_top_km_bubbles(bubbles, trim):
	quantile = get_km_quantile(bubbles,1-trim)
	
	trim_bubble_ids = []
	
	for bubble_id, bubble in bubbles.items():
		if bubble.allele0.km > quantile or bubble.allele1.km > quantile:
			trim_bubble_ids.append(bubble_id)
		
	for bubble_id in trim_bubble_ids:
		del bubbles[bubble_id]
		
		
def iterative_haplo_clust(bubbles, long_reads, q, itr=2, trim_km=0.05, min_haplotagged_bubbles = 1):
	
	'''
	Assign haplotypes to long reads and bubbles in {itr} iterations, starting from the first set of phased bubbles from the input file
	
	Parameters:
		itr (int)						: the number of halotype clustering iterations
		bubble_first_itr_phase_file: Path to the file containing the first set of phased bubbles
		bubbles			 				: A dictionary {bubble_id -> bubbles}
		long_reads		 				: A dictionary {long_read_name -> long_read}
	'''
	
	# trimming the set of bubbles with the highest km values
	trim_top_km_bubbles(bubbles, trim_km)
	
	#add_bubble_allele_pred_haplo(bubble_first_itr_phase_file, bubbles)
	add_bubbles_het_positions(bubbles)
	
	print('haploclust long reads first iteration...')
	print('setting alignments edit distance...')
	num = 0
	total_reads = len(long_reads)-(len(long_reads)%10)
	start_time = time.time()
	for read_name, long_read in long_reads.items():
		num += 1
		if num*10 % total_reads == 0:
			print(num*100/total_reads, '% of reads processed')
			
		if itr == 1:
			min_haplotagged_bubbles = 1
		
		long_read.link_alt_alignments()
		long_read.set_alignments_edit_dist(q)
		long_read.phase(min_haplotagged_bubbles = min_haplotagged_bubbles)
		
		# testing:
		#if long_read.pred_haplo != None:
		#	for aln in long_read.alignments:
		#		bubble_allele = aln.bubble_allele
		#		if bubble_allele.pred_haplo == None and bubble_allele.id == 0:
		#			test_bubbles[bubble_allele.bubble.id] = bubble_allele.bubble
					
	print('elapsed time =', time.time()-start_time)
	
	for it in range(itr-1):
		#TODO: the bubble phase eval here is different from first_itr dir output. They should be the same!!!
		print('outputting haploclust iteration', it+1)
		
		if it == itr-1:
			min_haplotagged_bubbles = 1
			
		print('haploclust iteration', it+2)
		haploclust_bubbles(bubbles)		
		haploclust_long_reads(long_reads, min_haplotagged_bubbles)


def output_phasing(bubbles, long_reads, bubble_phasing_file, long_read_phasing_file):
	print('outputting bubbles phasing')
	with open(bubble_phasing_file, 'w') as out:
		print("bubble_id\tallele0haplo", file=out)
		for bubble_id, bubble in bubbles.items():
			print(bubble_id, "\t", bubble.allele0.pred_haplo, file=out)

	print('outputting long reads phasing')
	with open(long_read_phasing_file, 'w') as out:
		print("long_read_name\thaplotype", file=out)
		for read_name, long_read in long_reads.items():
			print(read_name, "\t", long_read.pred_haplo, file=out)


