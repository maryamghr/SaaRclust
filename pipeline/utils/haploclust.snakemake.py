from haploclust import *
from bubble_long_read_alignment import *
from parsing import *
from evaluate_haploclust import *
from argparse import ArgumentParser


if __name__ == "__main__":
	parser = ArgumentParser(description=__doc__)
	parser.add_argument("--bubble_fasta_file", type=str, help="Bubble fasta file", required=True)
	parser.add_argument("--minimap_file", type=str, help="The set of minimap files containing alignments of bubbles to long reads", required=True)
	parser.add_argument("--long_reads_fasta_file", type=str, help="The set of long reads fasta files", required=True)
	parser.add_argument("--bubble_first_itr_phase_file", type=str, help="Bubble first iteration phase file", required=True)
	parser.add_argument("--bubble_phase_file", type=str, help="output bubble phase file", required=True)
	parser.add_argument("--long_read_phase_file", type=str, help="output long reads phase file", required=True)
	parser.add_argument("--itr", type=int, help="number of iterations for haplotype clustering", required=True)
	parser.add_argument("--min_bubbles", type=int, help="minimum number of heterozygous bubbles for haplotype clustering", required=True)
	parser.add_argument("--het_kmer_len", type=int, help="The length of heterozygous kmer for computing edit distance", required=True)
	parser.add_argument("--with_km", action='store_true', help="True if bubbles have km information in their name")

	
	args = parser.parse_args()

	with_km = True if "with_km" in args else False
	print('with_km', 	with_km)
	
	bubbles = get_bubbles(args.bubble_fasta_file, with_km)
	long_reads = get_long_reads(args.long_reads_fasta_file)
	set_alignments_from_minimap_file(args.minimap_file, bubbles, long_reads)
	
	add_bubble_allele_pred_haplo(args.bubble_first_itr_phase_file, bubbles)
	
	iterative_haplo_clust(bubbles=bubbles, long_reads=long_reads, q=args.het_kmer_len, itr=args.itr, min_haplotagged_bubbles=args.min_bubbles)

	output_phasing(bubbles, long_reads, args.bubble_phase_file, args.long_read_phase_file)
