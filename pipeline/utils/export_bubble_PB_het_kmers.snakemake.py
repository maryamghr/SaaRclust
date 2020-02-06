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

output_bubble_and_pb_kmers(snakemake.input["bubble_PB_minimap"], bubble_het_positions, bubble_allele_to_kmers, pb_name_to_seq, q, snakemake.output[0])


