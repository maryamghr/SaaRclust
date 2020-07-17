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


def get_eval_file_names(bubbles_haploclust_evaluation_file, long_reads_haploclust_evaluation_file):	
	bubble_eval_file_name_split = bubbles_haploclust_evaluation_file.split("iteration")
	bubble_eval_file_prefix = bubble_eval_file_name_split[0]+"iteration"
	bubble_eval_file_suffix = bubble_eval_file_name_split[1].split("_")
	bubble_eval_file_suffix = "_" + "_".join(bubble_eval_file_suffix[1:])
	
	bubble_eval_file_name = {"prefix": bubble_eval_file_prefix, "suffix": bubble_eval_file_suffix}
	
	long_read_eval_file_name_split = long_reads_haploclust_evaluation_file.split("iteration")
	long_read_eval_file_prefix = long_read_eval_file_name_split[0]+"iteration"
	long_read_eval_file_suffix = long_read_eval_file_name_split[1].split("_")
	long_read_eval_file_suffix = "_" + "_".join(long_read_eval_file_suffix[1:])
	
	long_read_eval_file_name = {"prefix": long_read_eval_file_prefix, "suffix": long_read_eval_file_suffix}
	
	return bubble_eval_file_name, long_read_eval_file_name
	
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
	parser.add_argument("--long_read_clust_files", nargs='*', help="The set of long reads chromosome clustered files", required=True)	
	parser.add_argument("--clust_to_chrom_file", type=str, help="Cluster to chrom mapping file", required=True)
	parser.add_argument("--long_read_haplotagged_bam_files", nargs='*', help="The set of long reads haplotagged bam files", required=True)	
	parser.add_argument("--bubble_phase_file", type=str, help="output bubble phase file", required=True)
	parser.add_argument("--test_bubble_phase_file", type=str, help="output test bubble phase file", required=True)
	parser.add_argument("--long_read_phase_file", type=str, help="output long reads phase file", required=True)
	parser.add_argument("--bubbles_haploclust_evaluation_file", type=str, help="The output bubbles clustring evaluation file", required=True)
	parser.add_argument("--bubbles_first_itr_haploclust_evaluation_file", type=str, help="The output first iteration bubbles clustring evaluation file", required=True)
	parser.add_argument("--bubbles_haplo_edit_dist_file", type=str, help="The output bubbles haplo edit dist file", required=True)
	parser.add_argument("--long_reads_haploclust_evaluation_file", type=str, help="The output long reads phasing evaluation file", required=True)
	parser.add_argument("--long_reads_haplo_edit_dist_file", type=str, help="The output long reads haplo edit dist file", required=True)
	parser.add_argument("--kmers_file", type=str, help="The output kmers file", required=True)
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
	
	add_bubble_allele_pred_haplo(args.bubble_first_itr_phase_file, bubbles)
	
	#if evaluation mode:
	add_bubbles_true_info(args.bubble_haplotagged_bam_file, bubbles)
	clust_to_chrom = get_clust_to_chrom(args.clust_to_chrom_file)
	add_bubble_clust(args.bubble_clust_file, bubbles)
	add_long_reads_true_info(args.long_read_haplotagged_bam_files, long_reads)
	
	evaluate_bubble_clustering(bubbles, clust_to_chrom, args.bubbles_first_itr_haploclust_evaluation_file)
	
	bubble_eval_file_name, long_read_eval_file_name = get_eval_file_names(args.bubbles_haploclust_evaluation_file, args.long_reads_haploclust_evaluation_file)

	print(bubble_eval_file_name, long_read_eval_file_name)	
	
	test_bubbles = iterative_haplo_clust(bubbles, long_reads, q, clust_to_chrom, bubble_eval_file_name, long_read_eval_file_name, itr)
	
	output_phasing(bubbles, long_reads, args.bubble_phase_file, args.long_read_phase_file)
	output_test_bubbles(test_bubbles, args.test_bubble_phase_file)

#if evaluation mode:
	evaluate_bubble_clustering(bubbles, clust_to_chrom, args.bubbles_haploclust_evaluation_file)
	output_bubbles_haplo_dist(bubbles, args.bubbles_haplo_edit_dist_file, with_km)
	evaluate_long_read_clustering(long_reads, args.long_reads_haploclust_evaluation_file)
	output_long_reads_haplo_dist(long_reads, args.long_reads_haplo_edit_dist_file)
	# FIXME: It takes a super long time after 100% !!!
	output_kmers(long_reads, args.kmers_file)
	
	#output_sampled_long_reads(100, (0.35, 0.4), args.long_reads_with_peak_frac_haplo_edit_dist)
	#output_sampled_long_reads(10, (0, 0.1), args.long_reads_with_small_frac_haplo_edit_dist)
	
	