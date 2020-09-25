#from output_valid_maps import *
#import pdb
from parsing import *

#libs = snakemake.params['libs']
#
#libs_true_order = True
#for i in range(len(libs)):
#	print(libs[i], snakemake.input['ss_reads'][i])
#	print(libs[i], snakemake.input['map'][i])
#
#	if snakemake.input['ss_reads'][i].find(libs[i])==-1 or snakemake.input['map'][i].find(libs[i])==-1:
#		libs_true_order = False
#
#print('libs_true_order:', libs_true_order)

#if snakemake.params['input_type']=='unitig':
#	unitig_to_bubble_allele = map_unitig_to_bubble_allele(snakemake.input['bubbles'])

#clust_pair = snakemake.wildcards["clust_pair"].split('_')

#output_valid_maps(snakemake.input['ss_reads'], snakemake.input['ss_clust_file'], clust_pair, snakemake.input['unitigs'], snakemake.input['map'], snakemake.output[0], snakemake.log[0], snakemake.params['libs'], snakemake.params['input_type'], unitig_to_bubble_allele)

unitig_to_bubble_allele = map_unitig_to_bubble_allele(snakemake.input['bubbles'])
output_bwa_fastmap_matches(unitig_to_bubble_allele, snakemake.input['map'], snakemake.output[0])