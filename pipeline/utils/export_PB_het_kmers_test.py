import unittest
from export_PB_het_kmers import *


class Testexporthetkmers(unittest.TestCase):
	
	def test_reverse_comp(self):
		self.assertEqual('acGT', reversecomp('ACgt'))
		
	def test_reverse_comp_empty_input(self):
		self.assertEqual('', reversecomp(''))
	
	def test_add_bubble_kmers(self):
		bubble_het_positions, bubble_allele_to_kmers = {}, {}
		add_bubble_kmers(bubble_het_positions, bubble_allele_to_kmers, 12, 'AGCTC', 'AGGTC', 1)
		self.assertEqual(bubble_het_positions, {12:[2]})
		self.assertEqual(bubble_allele_to_kmers, {(12,0):['GCT'], (12,1):['GGT']}) 
		
	def test_add_kmers_bubbles_with_two_het_pos(self):
		bubble_het_positions, bubble_allele_to_kmers = {}, {}
		add_bubble_kmers(bubble_het_positions, bubble_allele_to_kmers, 12, 'AGCTCA', 'AGGTTA', 1)
		self.assertEqual(bubble_het_positions, {12:[2,4]})
		self.assertEqual(bubble_allele_to_kmers, {(12,0):['GCT','TCA'], (12,1):['GGT','TTA']})


	def test_find_query_interval(self):
		ref_start, ref_end, alignment_start_pos = 4, 5, 2
		cigar='2M1D2M2I1X'
		aligned_ref_positions   = [2,3, 4, 5,6,'-','-',7]
		aligned_query_positions = [0,1,'-',2,3, 4 , 5 ,6]
		self.assertEqual(find_query_interval(ref_start, ref_end, alignment_start_pos, cigar), (2,2))
		
	def test_find_reference_interval(self):
		query_start, query_end, aln_ref_start_pos, aln_query_start_pos = 2, 5, 2, 0
		cigar='2M1D2M2I1X'
		aligned_ref_positions   = [2,3, 4, 5,6,'-','-',7]
		aligned_query_positions = [0,1,'-',2,3, 4 , 5 ,6]
		self.assertEqual(find_reference_interval(query_start, query_end, aln_ref_start_pos, aln_query_start_pos, cigar), (5,6))
		
	def test_read_het_snv_bubbles(self):
		bubbles_file, q = 'example_files/example_bubbles_withclust.fa', 10
		bubble_het_positions, bubble_allele_to_kmers = read_het_snv_bubbles(bubbles_file, q)
		self.assertEqual(bubble_het_positions, {2:[62]})
		self.assertEqual(bubble_allele_to_kmers, {(2,0):['ATAAATGCCAGGAATTCTGTA'], (2,1):['ATAAATGCCAAGAATTCTGTA']})
		
	def test_read_het_snv_bubbles_with_last_clustered_bubbles(self):
		bubbles_file, q = 'example_files/example_bubbles_with_clustered_last_bubble.fa', 10
		bubble_het_positions, bubble_allele_to_kmers = read_het_snv_bubbles(bubbles_file, q)
		self.assertEqual(bubble_het_positions, {2:[62], 3:[62]})
		self.assertEqual(bubble_allele_to_kmers, {(2,0):['ATAAATGCCAGGAATTCTGTA'], (2,1):['ATAAATGCCAAGAATTCTGTA'], 
			(3,0):['AAAAGAGAGGAACAAATGGCA'], (3,1):['AAAAGAGAGGTACAAATGGCA']})
	
