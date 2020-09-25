from bubble_long_read_alignment import *
import pysam
import time
import pdb
import gzip
from argparse import ArgumentParser

def get_alignments_from_minimap_file(minimap_files_list):

	'''
	Given a list of minimap alignment files, and bubbles and long_reads, creates alignment objects
	
	Parameters:
		minimap_files_list: A list of paths to minimap alignment files
		bubbles			 : A dictionary {bubble_id -> bubbles}
		long_reads		 : A dictionary {long_read_name -> long_read}
	'''
	
	bubbles = {}
	long_reads = {}

	for minimap_file in minimap_files_list:
		start_time = time.time()
		print('reading alignments from file', minimap_file)
		with gzip.open(minimap_file) as minimap:
			
			for line in minimap:
				line = line.decode("utf-8")
				sp = line.split()
					
				bubble_name, bubble_len, bubble_start, bubble_end, strand, \
				read_name, long_read_len, long_read_start, long_read_end = \
				sp[0], int(sp[1]), int(sp[2]), int(sp[3])-1, sp[4], \
				sp[5], int(sp[6]), int(sp[7]), int(sp[8])-1
				
				bubble_name_sp = bubble_name.split('_')
				
				bubble_id, bubble_allele_id = int(bubble_name_sp[1]), int(bubble_name_sp[3])-1
				
				assert(bubble_allele_id == 0 or bubble_allele_id == 1), 'bubble ' + str(bubble_id) + ': allele should be 0 or 1'

				# remove the beginning part of the cigar string
				cigar = sp[-1] # assuming that always the last tag is cigar
				cigar = cigar.split('cg:Z:')[1]
				
				read_name_sp = read_name.split('/ccs')
				read_name = read_name_sp[0]+'/ccs'

				if bubble_id in bubbles:
					bubble = bubbles[bubble_id]
				else:
					bubble = Bubble(bubble_id)
					bubbles[bubble_id] = bubble
					
				bubble_allele = bubble.allele0 if bubble_allele_id == 0 else bubble.allele1
				
				if bubble_allele == None:
					bubble_allele = BubbleAllele(bubble_allele_id, bubble_name)
					bubble.add_allele(bubble_allele)
					
				if read_name in long_reads:		
					long_read = long_reads[read_name]
				else:
					long_read = LongRead(read_name)
					long_reads[read_name] = long_read				
				
				aln = Alignment(long_read=long_read, bubble_allele=bubble_allele)
				
		print('elapsed time =', time.time()-start_time)

	pdb.set_trace()


if __name__ == "__main__":
	parser = ArgumentParser(description=__doc__)
	parser.add_argument("--minimap_files_list", nargs='*', help="The set of minimap files containing alignments of bubbles to long reads", required=True)
	
	args = parser.parse_args()

	get_alignments_from_minimap_file(args.minimap_files_list)