from evaluate_bubbles_phase import *

num_bubbles, bubble_id_to_chrom, ground_true_bubble_id_to_h0_allele = get_ground_true_bubbles_chrom_and_h0allele(snakemake.input["bubble_haplotagged_bam_file"])
clust_to_chrom = get_clust_to_chrom(snakemake.input["clust_to_chrom_file"])
bubble_id_to_clust = get_bubble_id_to_clust(snakemake.input["bubble_clust_file"])
bubble_id_to_h0_allele = get_bubble_id_to_h0_allele(snakemake.input["bubble_phase_file"])

print_dict_head(bubble_id_to_chrom, 10)
print_dict_head(bubble_id_to_h0_allele, 10)
print_dict_head(bubble_id_to_clust, 10)

print('len(bubble_id_to_chrom) =', len(bubble_id_to_chrom), ', len(bubble_id_to_h0_allele)=', len(bubble_id_to_h0_allele), ', len(bubble_id_to_clust)=', len(bubble_id_to_clust))

evaluate_bubble_phase(num_bubbles, bubble_id_to_chrom, ground_true_bubble_id_to_h0_allele, clust_to_chrom, bubble_id_to_clust, bubble_id_to_h0_allele, snakemake.output["evaluation"], snakemake.output["false_phased_bubbles_file"])

