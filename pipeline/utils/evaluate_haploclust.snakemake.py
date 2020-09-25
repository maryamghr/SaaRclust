from haploclust import *
from bubble_long_read_alignment import *
from parsing import *
from evaluate_haploclust import *
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("--bubble_haplotagged_bam_file", type=str, help="Bubble haplotagged bam file", required=True)
parser.add_argument("--bubble_clust_file", type=str, help="Bubbles cluster file", required=True)
parser.add_argument("--clust_to_chrom_file", type=str, help="Cluster to chrom mapping file", required=True)
parser.add_argument("--long_read_haplotagged_bam_files", nargs='*', help="The set of long reads haplotagged bam files", required=True)	
parser.add_argument("--bubble_phase_files", nargs='*', help="output bubble phase file", required=True)
parser.add_argument("--long_read_phase_files", nargs='*', help="output long reads phase file", required=True)
parser.add_argument("--long_reads_clust_files", nargs='*', help="output long reads phase file", required=True)
parser.add_argument("--bubbles_haploclust_evaluation_file", type=str, help="The output bubbles clustring evaluation file", required=True)
parser.add_argument("--long_reads_haploclust_evaluation_file", type=str, help="The output long reads phasing evaluation file", required=True)
	
args = parser.parse_args()

bubbles = get_bubbles_from_bam(args.bubble_haplotagged_bam_file)
long_reads = get_long_reads_from_bam(args.long_read_haplotagged_bam_files)
	
for bubble_phase_file in args.bubble_phase_files:
	add_bubble_allele_pred_haplo(bubble_phase_file, bubbles)
add_long_reads_pred_haplotype(args.long_read_phase_files, long_reads)
	
clust_to_chrom = get_clust_to_chrom(args.clust_to_chrom_file)
add_bubble_clust(args.bubble_clust_file, bubbles)
add_long_reads_clust(args.long_reads_clust_files, long_reads)
	
evaluate_bubble_clustering(bubbles, clust_to_chrom, args.bubbles_haploclust_evaluation_file)
evaluate_long_read_clustering(long_reads, args.long_reads_haploclust_evaluation_file)
	

