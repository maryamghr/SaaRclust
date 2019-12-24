from evaluate_long_reads_phasing import *

print("snakemake.input[\"long_read_haplotagged_bam_files\"]")
print(snakemake.input["long_read_haplotagged_bam_files"])
print("type of input")
print(type(snakemake.input["long_read_haplotagged_bam_files"]))
chrom_to_num_ground_true_phased_reads, ground_true_read_to_chrom, ground_true_read_to_haplo = get_ground_true_chrom_haplo(snakemake.input["long_read_haplotagged_bam_files"])
read_to_haplo = get_reads_haplotypes(snakemake.input["long_read_phase_file"])
output_long_reads_phasing_accuracy(chrom_to_num_ground_true_phased_reads, ground_true_read_to_chrom, ground_true_read_to_haplo,  read_to_haplo, snakemake.output["evaluation"])

