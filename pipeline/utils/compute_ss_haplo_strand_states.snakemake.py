from compute_ss_haplo_strand_states import *

chrom_dir_to_clust = get_chrom_dir_to_clust_map(snakemake.input["clust_partners_file"])
ss_counts_in_chroms = get_ss_counts_in_chroms(snakemake.input["ss_counts_file"])
output_ss_haplo_strand_states(chrom_dir_to_clust, ss_counts_in_chroms, snakemake.output["SS_haplo_strand_states"], snakemake.params["min_w_frac_in_wc_state"], snakemake.params["max_w_frac_in_wc_state"], snakemake.params["max_haplo_count_ratio"])
