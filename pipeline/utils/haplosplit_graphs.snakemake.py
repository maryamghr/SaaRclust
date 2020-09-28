from process_graph import *

#haplo_file = '../aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/split/hifiasm/phased_bubbles/iteration1_V25_V32_r.data'
#input_gfa = '../aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/split/hifiasm/V25_V32.r_utg.gfa'
#h1_gfa = '../aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/split/hifiasm/V25_V32.r_h1_utg.gfa'
#h2_gfa = '../aligns_k15_w1_f0.1_z500/SaaRclust_results_HG00733/Clusters/split/hifiasm/V25_V32.r_h2_utg.gfa'

node_to_haplo = get_node_haplotypes(snakemake.input['haplo_file'])
read_haplo = separate_haplotypes(node_to_haplo, snakemake.input['input_gfa'], snakemake.output['h1_gfa'], snakemake.output['h2_gfa'])
output_read_haplo(read_haplo, snakemake.input['ccs_fasta'], snakemake.output['ccs_phase'])