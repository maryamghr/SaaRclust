class Bubble:

	def __init__(self, id=None):
		self.id = id
		self.actual_chrom = None
		self.clust = None
		self.allele0, self.allele1 = None, None
		self.actual_type = None # in ["unmapped", "invalid_chrom", "untagged"]
		self.pred_type = None # in ["true_haplo", "false_haplo", "haploclust_false_pos", "not_chrom_clust", "garbage_clust", "not_haplo_clust"]

	def print(self):
		print('id =', self.id)
		print('actual_chrom =', self.actual_chrom)
		print('clust =', self.clust)
		print('actual_type =', self.actual_type)
		print('pred_type =', self.pred_type)
		print('****************************\nallele0\n****************************')
		if self.allele0 != None:
			self.allele0.print()
		else:
			print('allele0: None')

		print('****************************\nallele1\n****************************')
		if self.allele1 != None:
			self.allele1.print()
		else:
			print('allele1: None')

		print('****************************')

	def add_allele(self, bubble_allele):
		assert (bubble_allele.id==0 or bubble_allele.id==1), 'bubble allele =, ' + str(bubble_allele.id) + 'should be 0 or 1'
		if bubble_allele.id == 0:
			self.allele0 = bubble_allele
		else:
			self.allele1 = bubble_allele

	def get_haplo_num_aligned_reads(self):
		haplo_num_aln_reads = (None, None)

		if self.allele0.pred_haplo != None:
			haplo_num_aln_reads = (len(self.allele0.alignments), len(self.allele1.alignments)) if self.allele0.pred_haplo==0 else \
			       								(len(self.allele1.alignments), len(self.allele0.alignments))

		return haplo_num_aln_reads



class BubbleAllele:

	def __init__(self, id=None, bubble=None):
		self.id = id
		self.bubble = bubble
		self.actual_haplo, self.pred_haplo = None, None
		#self.haplo_edit_dist = [0,0]
		self.km = 0
		self.alignments = []

	def print(self):
		print('allele id =', self.id)
		print('bubble_id =', self.bubble.id)
		print('km =', self.km)
		print('actual_haplo =', self.actual_haplo)
		print('pred_haplo =', self.pred_haplo)
		print('alignments =', self.alignments)

	def get_haplotypes_edit_dist(self):

		haplo_edit_dist = [0,0]

		for aln in self.alignments:
			long_read_haplo = aln.long_read.pred_haplo

			if long_read_haplo == None:
				continue

			assert (long_read_haplo==0 or long_read_haplo==1), 'the predicted haplotype for long read ' + \
				aln.long_read.name + ' is ' + str(long_read_haplo) + ', should be 0 or 1'

			haplo_edit_dist[long_read_haplo] += aln.edit_dist

		return haplo_edit_dist


class LongRead:

	def __init__(self, name):
		self.name = name
		self.actual_haplo, self.pred_haplo = None, None
		self.actual_chrom, self.clust = None, None
		self.haplo0_edit_dist, self.haplo1_edit_dist = 0, 0
		self.actual_type, self.pred_type = None, None
		self.alignments = []

	def set_haplotypes_edit_dist(self):

		haplo_edit_dist = [0, 0] # h0_dist, h1_dist, respectively

		for aln in self.alignments:
			bubble_allele_haplo = aln.bubble_allele.pred_haplo

			if bubble_allele_haplo == None:
				continue

			assert (bubble_allele_haplo==0 or bubble_allele_haplo==1), 'error in alignment ' + self.name + \
				', bubble ' + str(aln.bubble_allele.bubble.id) + ', allele ' + \
				str(aln.bubble_allele.id) + ': bubble allele predicted haplotype is ' + \
				str(bubble_allele_haplo) + ', should be 0 or 1'

			haplo_edit_dist[bubble_allele_haplo] += aln.edit_dist

		self.haplo0_edit_dist, self.haplo1_edit_dist = haplo_edit_dist[0], haplo_edit_dist[1]
		

class Alignment:
	
	def __init__(self, long_read, bubble_allele, edit_dist):
		self.long_read, self.bubble_allele, self.edit_dist = long_read, bubble_allele, edit_dist
		self.long_read.alignments.append(self)
		self.bubble_allele.alignments.append(self)

	def __init__(self, long_read, bubble_allele, bubble_kmer, long_read_kmer, edit_dist):
		self.long_read, self.bubble_allele = long_read, bubble_allele
		self.long_read.alignments.append(self)
		self.bubble_allele.alignments.append(self)
		self.bubble_kmer = bubble_kmer
		self.long_read_kmer = long_read_kmer
		self.edit_dist = edit_dist
		
	def output_print(self):
		# bubbleName	bubbleAllele	PBname	bubbleKmer	PBkmer	kmersEditDistance bubble_alle_pred_haplo	long_read_pred_haplo
		print_str = str(self.bubble_allele.bubble.id)
		print_str += '\t' + str(self.bubble_allele.id)
		print_str += '\t' + str(self.long_read.name)
		print_str += '\t' + str(self.bubble_kmer)
		print_str += '\t' + str(self.long_read_kmer)
		print_str += '\t' + str(self.edit_dist)
		print_str += '\t' + str(self.bubble_allele.pred_haplo)
		print_str += '\t' + str(self.long_read.pred_haplo)
		
		return print_str
		
