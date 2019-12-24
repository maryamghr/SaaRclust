from count_ss_reads import *

output_ss_counts_in_chroms(snakemake.input["ss_haplotagged_bam"], snakemake.output["ss_counts_file"])
