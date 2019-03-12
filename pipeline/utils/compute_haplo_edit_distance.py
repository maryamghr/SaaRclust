from whatshap.align import edit_distance
import sys

with open(snakemake.input[0]) as f:
	with open(snakemake.output[0], 'w') as out:
		for line in f:
			if line=="":
				break
		
			sp = line.split()

			if sp[0] == "PB_name":
				print(line.strip() + "\th1_dist\th2_dist", file=out)
				continue

			pb_name = sp[0]
			pb_seq = sp[1]
			h1_kmer = sp[2]
			h2_kmer = sp[3]
			
			d1 = edit_distance(pb_seq, h1_kmer)
			d2 = edit_distance(pb_seq, h2_kmer)

			print(line.strip()+"\t"+str(d1)+"\t"+str(d2), file=out)
