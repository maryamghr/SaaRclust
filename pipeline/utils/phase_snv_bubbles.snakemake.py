from phase_snv_bubbles import *

print("getting haplo strand states ...")
lib_clust_to_haplo = read_strand_states(list(snakemake.input["SS_haplo_strand_states"]))

print("phasing the bubbles and writing the phase information in the output file ...")
phase_bubbles(snakemake.input["SS_bubble_map"], lib_clust_to_haplo, snakemake.output[0])
