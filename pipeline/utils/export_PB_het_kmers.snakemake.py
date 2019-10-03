import sys
import copy
import gzip
from export_PB_het_kmers import *
####################################
# Note: q should be shorter than the kmer length...
q = snakemake.params["het_kmer_len"]
print('q =', q)
bubble_het_positions, bubble_allele_to_kmers = read_het_snv_bubbles(snakemake.input["bubbles"], q)
print('len(bubble_het_positions) =', len(bubble_het_positions), ', len(bubble_het_positions) =', len(bubble_het_positions))
print('head(bubble_het_positions) =', print_dict_head(bubble_het_positions, 5))
print('head(bubble_het_positions) =', print_dict_head(bubble_het_positions, 5))
ss_to_kmer_altkmer, ss_to_kmer_pos_and_interval, ss_to_clust = read_SS_bubble_map(list(snakemake.input["SS_bubble_map"]), bubble_het_positions, bubble_allele_to_kmers, q)
print('len(ss_to_kmer_altkmer) =', len(ss_to_kmer_altkmer), ', len(ss_to_kmer_pos_and_interval) =', len(ss_to_kmer_pos_and_interval), ', len(ss_to_clust) =', len(ss_to_clust))
print('head(ss_to_kmer_altkmer) =', print_dict_head(ss_to_kmer_altkmer, 5))
print('head(ss_to_kmer_pos_and_interval) =', print_dict_head(ss_to_kmer_pos_and_interval, 5))
print('head(ss_to_clust) =', print_dict_head(ss_to_clust, 5))
lib_clust_to_haplo = read_strand_states(list(snakemake.input["SS_haplo_strand_states"]))
print('len(lib_clust_to_haplo) =', len(lib_clust_to_haplo))
print('head(lib_clust_to_haplo) =', print_dict_head(lib_clust_to_haplo, 5))
pb_to_kmers, pb_to_kmer_intervals = get_pb_kmers(snakemake.input["SS_PB_minimap"], ss_to_kmer_altkmer, ss_to_kmer_pos_and_interval, ss_to_clust, lib_clust_to_haplo)
print('len(pb_to_kmers) =', len(pb_to_kmers), 'len(pb_to_kmer_intervalsprint) =', len(pb_to_kmer_intervals))
output_pb_het_kmers(snakemake.input["PB_fasta"], snakemake.output[0], pb_to_kmers, pb_to_kmer_intervals)
