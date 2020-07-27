from export_bubble_PB_het_kmers import *
####################################
# Note: q should be shorter than the kmer length...

q = snakemake.params["het_kmer_len"]
print('q =', q)

bubble_het_positions, bubble_allele_to_kmers = read_het_snv_bubbles(snakemake.input["bubbles"], snakemake.input["bubbles_clust_file"], q)

print('len(bubble_het_positions) =', len(bubble_het_positions), ', len(bubble_allele_to_kmers) =', len(bubble_allele_to_kmers))
print('head(bubble_het_positions) =')
print_dict_head(bubble_het_positions, 5)
print('head(bubble_allele_to_kmers) =')
print_dict_head(bubble_allele_to_kmers, 5)

pb_name_to_seq = get_pb_name_to_seq(snakemake.input["PB_fasta"])

print('head(pb_name_to_seq) =')
print_dict_head(pb_name_to_seq, 5)

output_bubble_and_pb_kmers(snakemake.input["bubble_PB_minimap"], snakemake.params["bubble_info"], bubble_het_positions, bubble_allele_to_kmers, pb_name_to_seq, q, snakemake.output[0])


# run in cmd:
#q=10                                                                                                       
#bubbles_file, bubbles_clust_file = '../bubbles/snv_bubbles_k63_a3_l23.fa', '../bubbles/snv_bubbles_clusters_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00733.data'
#PB_fasta, bubble_PB_minimap = '../../../ccs_bams/HG00733.000.fasta', '../aligns_k15_w1_f0.1_z500/snv_bubbles_k63_a3_l23/HG00733_000.paf.gz'
#bubble_het_positions, bubble_allele_to_kmers = read_het_snv_bubbles(bubbles_file, bubbles_clust_file, q)

