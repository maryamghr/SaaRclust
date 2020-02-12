from split_phased_long_reads_by_cluster import *

split_phase_files_by_clust(snakemake.input['phase_files'], snakemake.input['clust_files'], snakemake.output)
