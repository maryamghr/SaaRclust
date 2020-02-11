from phase_long_reads import *

update_bubble_phase(snakemake.input['het_kmers_edit_dist'], snakemake.input['long_reads_phase_file'], snakemake.output[0])
