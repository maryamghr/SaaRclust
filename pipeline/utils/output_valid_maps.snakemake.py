from output_valid_maps import *
import pdb

if snakemake.params['input_type']=='unitig':
	unitig_to_bubble_allele = map_unitig_to_bubble_allele(snakemake.input['bubbles'])

output_valid_maps(snakemake.input['ss_reads'], snakemake.input['unitigs'], snakemake.input['map'], snakemake.output[0], snakemake.log[0], snakemake.wildcards['lib'], snakemake.params['input_type'], unitig_to_bubble_allele)
