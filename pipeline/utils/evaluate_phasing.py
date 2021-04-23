#from haploclust import *
from bubble_long_read_alignment import *
from parsing import *
from evaluate_haploclust import *
from argparse import ArgumentParser
import pdb

parser = ArgumentParser(description=__doc__)
#parser.add_argument("-b", help="The type of the sequence is bubble")
parser.add_argument("--haplotagged_bam_files", nargs="*", help="haplotagged bam file", required=True)
parser.add_argument("--phase_files", nargs="*", help="phase file", required=True)
parser.add_argument("--output_file", type=str, help="output file", required=True)
	
args = parser.parse_args()

#bubbles = get_bubbles_from_bam(args.haplotagged_bam_file)
#print_reference_mapping_stats(bubbles)
#if 'b' in args:
#add_bubble_allele_pred_haplo(args.phase_file, bubbles)
#evaluate_phasing(bubbles)

long_reads = get_long_reads_from_bam(args.haplotagged_bam_files)
add_long_reads_pred_haplotype(args.phase_files, long_reads)
pdb.set_trace()
evaluate_long_read_clustering(long_reads, args.output_file)
