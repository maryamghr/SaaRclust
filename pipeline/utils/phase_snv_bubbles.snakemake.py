from phase_snv_bubbles import *

print("getting haplo strand states ...")
lib_clust_to_haplo = read_strandphaser_strand_states(snakemake.input["phased_strand_states"])

print("getting ss clusters ...")
ss_to_clust = get_ss_clust(snakemake.input["ss_clust"])

print("phasing the bubbles and writing the phase information in the output file ...")
phase_bubbles(snakemake.input["SS_bubble_map"], ss_to_clust, lib_clust_to_haplo, snakemake.output[0])
