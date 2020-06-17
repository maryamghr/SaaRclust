from haploclust import *
from bubble_long_read_alignment import *
from parsing import *
from evaluate_haploclust import *
from argparse import ArgumentParser

#############################################################################################
#with_km=snakemake.params['with_km']
#q=int(snakemake.params['het_kmer_len'])
#itr=int(snakemake.params['itr'])
	
#bubbles = get_bubbles(snakemake.input["bubble_fasta_file"], with_km)
#long_reads = get_long_reads(snakemake.input["long_reads_fasta_files"])
#set_alignments_from_minimap_file(snakemake.input["minimap_files_list"], bubbles, long_reads)
	
#iterative_haplo_clust(itr, snakemake.input["bubble_first_itr_phase_file"], bubbles, long_reads, q)
	
#output_phasing(snakemake.input["bubbles_phase_file"], snakemake.input["long_reads_phase_file"])
	
#if snakemake.params["evaluation"]:
#	add_bubbles_true_info(snakemake.input["bubble_haplotagged_bam_file"], bubbles)
#	clust_to_chrom = get_clust_to_chrom(snakemake.input["clust_to_chrom_file"])
#	add_bubble_clust(snakemake.input["bubble_clust_file"], bubbles)
#	add_long_reads_true_info(snakemake.input["long_reads_haplotagged_bam_files"], long_reads)
#	evaluate_bubble_clustering(bubbles, clust_to_chrom, snakemake.input["bubbles_clusering_evaluation_file"])
#	output_bubbles_haplo_dist(bubbles, snakemake.input["bubbles_haplo_edit_dist_file"], with_km)
#	evaluate_long_read_clustering(long_reads, snakemake.input["long_read_phase_evaluation_file"])
#	output_long_reads_haplo_dist(long_reads, snakemake.input["long_reads_haplo_edit_dist_file"])	
	#output_sampled_long_reads(100, (0.35, 0.4), args.long_reads_with_peak_frac_haplo_edit_dist)
	#output_sampled_long_reads(10, (0, 0.1), args.long_reads_with_small_frac_haplo_edit_dist)



with_km=True
q=10
itr=2

		
if __name__ == "__main__":
	parser = ArgumentParser(description=__doc__)
	parser.add_argument("--bubble_fasta_file", type=str, help="Bubble fasta file", required=True)
	parser.add_argument("--minimap_files_list", nargs='*', help="The set of minimap files containing alignments of bubbles to long reads", required=True)
	parser.add_argument("--long_reads_fasta_files", nargs='*', help="The set of long reads fasta files", required=True)
	parser.add_argument("--bubble_haplotagged_bam_file", type=str, help="Bubble haplotagged bam file", required=True)
	parser.add_argument("--bubble_first_itr_phase_file", type=str, help="Bubble first iteration phase file", required=True)
	parser.add_argument("--bubble_clust_file", type=str, help="Bubbles cluster file", required=True)
	parser.add_argument("--clust_to_chrom_file", type=str, help="Cluster to chrom mapping file", required=True)
	parser.add_argument("--long_read_haplotagged_bam_files", nargs='*', help="The set of long reads haplotagged bam files", required=True)	
	parser.add_argument("--bubble_phase_file", type=str, help="output bubble phase file", required=True)
	parser.add_argument("--long_read_phase_file", type=str, help="output long reads phase file", required=True)
	parser.add_argument("--bubbles_haploclust_evaluation_file", type=str, help="The output bubbles clustring evaluation file", required=True)
	parser.add_argument("--bubbles_haplo_edit_dist_file", type=str, help="The output bubbles haplo edit dist file", required=True)
	parser.add_argument("--long_reads_haploclust_evaluation_file", type=str, help="The output long reads phasing evaluation file", required=True)
	parser.add_argument("--long_reads_haplo_edit_dist_file", type=str, help="The output long reads haplo edit dist file", required=True)
	parser.add_argument("--itr", type=str, help="number of iterations for haplotype clustering", required=True)
	parser.add_argument("--het_kmer_len", type=str, help="The length of heterozygous kmer for computing edit distance", required=True)
	parser.add_argument("--with_km", action='store_true', help="True if bubbles have km information in their name")

	# testing output files for observation 
	#parser.add_argument("--long_reads_with_small_frac_haplo_edit_dist", type=str, help="The output sample long reads with small fraction of haplo_edit_dist", required=True)
	#parser.add_argument("--long_reads_with_peak_frac_haplo_edit_dist", type=str, help="The output sample long reads with peak fraction of haplo_edit_dist", required=True)
	######################################
	
	args = parser.parse_args()

	#counts(args.sample, args.input_bam, args.input_bed, args.counts_output)
	with_km = True if "with_km" in args else False
	print('with_km', 	with_km)
	
	q=int(args.het_kmer_len)
	itr=int(args.itr)

	bubbles = get_bubbles(args.bubble_fasta_file, with_km)
	long_reads = get_long_reads(args.long_reads_fasta_files)
	set_alignments_from_minimap_file(args.minimap_files_list, bubbles, long_reads)
	
	iterative_haplo_clust(args.bubble_first_itr_phase_file, bubbles, long_reads, q, itr)
	
	output_phasing(bubbles, long_reads, args.bubble_phase_file, args.long_read_phase_file)

#if evaluation mode:
	add_bubbles_true_info(args.bubble_haplotagged_bam_file, bubbles)
	clust_to_chrom = get_clust_to_chrom(args.clust_to_chrom_file)
	add_bubble_clust(args.bubble_clust_file, bubbles)
	add_long_reads_true_info(args.long_read_haplotagged_bam_files, long_reads)

	evaluate_bubble_clustering(bubbles, clust_to_chrom, args.bubbles_haploclust_evaluation_file)
	output_bubbles_haplo_dist(bubbles, args.bubbles_haplo_edit_dist_file, with_km)
	evaluate_long_read_clustering(long_reads, args.long_reads_haploclust_evaluation_file)
	output_long_reads_haplo_dist(long_reads, args.long_reads_haplo_edit_dist_file)	
	#output_sampled_long_reads(100, (0.35, 0.4), args.long_reads_with_peak_frac_haplo_edit_dist)
	#output_sampled_long_reads(10, (0, 0.1), args.long_reads_with_small_frac_haplo_edit_dist)
	
	