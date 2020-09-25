import matplotlib.pyplot as plt
from evaluate_long_reads_phasing import *

#def scatter_plot_h1_h2_dist(plot_file, long_reads_phase_file_list, het_reads_true_tagged_file, het_reads_false_tagged_file, het_reads_untagged_file):
	


# defining input files and calling the plotting function
chunkID = ['0'*(3-len(str(i)))+str(i) for i in range(100)]
long_reads_phase_file_list = ['../aligns_k15_w1_f0.1_z500/phased_long_reads/iteration2_pb_phase_chunk'+x+'_k63_a3_l23_HG00733.data' for x in chunkID]
het_reads_true_tagged_file='../evaluation/iteration2_long_read_phase_evaluation/het_reads_true_tagged_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00733.data'
het_reads_false_tagged_file='../evaluation/iteration2_long_read_phase_evaluation/het_reads_false_tagged_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00733.data'
het_reads_untagged_file='../evaluation/iteration2_long_read_phase_evaluation/het_reads_untagged_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00733.data'
plot_file='../evaluation/iteration2_long_read_phase_evaluation/plot_haplo_edit_dist_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00733.png'

read_to_haplo, read_to_haplo_dist = get_reads_haplotypes(long_reads_phase_file_list)
with open(het_reads_true_tagged_file) as f:
	true_tagged = [line.strip().split()[1] for line in f.readlines() if line!=""]
with open(het_reads_false_tagged_file) as f:
	false_tagged = [line.strip().split()[1] for line in f.readlines() if line!=""]
with open(het_reads_untagged_file) as f:
	untagged = [line.strip().split()[1] for line in f.readlines() if line!=""]

#true_tagged, false_tagged, untagged = [], [], []
#for chrom in valid_chroms:
#	true_tagged += het_reads_true_tagged[chrom]
#	false_tagged += het_reads_false_tagged[chrom]
#	untagged += het_reads_untagged[chrom]

true_tagged_haplo_dist = [read_to_haplo_dist[read] for read in true_tagged]
false_tagged_haplo_dist = [read_to_haplo_dist[read] for read in false_tagged]
untagged_haplo_dist = [read_to_haplo_dist[read] for read in untagged]

data = (true_tagged_haplo_dist, false_tagged_haplo_dist, untagged_haplo_dist)
colors = ("red", "green", "blue")
groups = ("true_tagged", "false_tagged", "untagged")

# Create plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, facecolor="1.0")

for data, color, group in zip(data, colors, groups):
	x, y = data
	ax.scatter(x, y, alpha=0.8, c=color, edgecolors='none', s=30, label=group)

plt.title('Matplot scatter plot')
plt.legend(loc=2)
plt.show()
plt.savefig(plot_file)
