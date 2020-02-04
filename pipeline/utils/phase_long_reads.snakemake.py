from phase_long_reads import *

phase_long_reads(snakemake.input['het_kmers_edit_dist'], snakemake.input['bubble_phase_file'], snakemake.output[0])
