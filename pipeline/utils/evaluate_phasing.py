from haploclust import *
from bubble_long_read_alignment import *
from parsing import *
from evaluate_haploclust import *
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
#parser.add_argument("-b", help="The type of the sequence is bubble")
parser.add_argument("--haplotagged_bam_file", type=str, help="Bubble haplotagged bam file", required=True)
parser.add_argument("--phase_file", type=str, help="bubble phase file", required=True)
	
args = parser.parse_args()

bubbles = get_bubbles_from_bam(args.haplotagged_bam_file)
print_reference_mapping_stats(bubbles)
#if 'b' in args:
add_bubble_allele_pred_haplo(args.phase_file, bubbles)
evaluate_phasing(bubbles)

#add_long_reads_pred_haplotype(args.long_read_phase_files, long_reads)
#evaluate_long_read_clustering(long_reads, args.long_reads_haploclust_evaluation_file)
	

