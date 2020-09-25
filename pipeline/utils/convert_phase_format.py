from __future__ import division
from bubble_long_read_alignment import *
from parsing import *
from argparse import ArgumentParser

if __name__ == "__main__":
	parser = ArgumentParser(description=__doc__)
	parser.add_argument("--bubble_fasta_file", type=str, help="Bubble fasta file", required=True)
	parser.add_argument("--input_phase_file", type=str, help="input bubble phase file", required=True)
	parser.add_argument("--output_phase_file", type=str, help="output bubble phase file", required=True)
	
	args = parser.parse_args()

	with_km = True
	bubbles = get_bubbles(args.bubble_fasta_file, with_km)
	
	with open(args.output_phase_file, 'w') as out:
		with open(args.input_phase_file) as f:
			
			next(f)
			
			for line in f:
				line_sp = line.split()
				
				bubble_id, pred_haplo = int(line_sp[0]), line_sp[1]				
				bubble = bubbles[bubble_id]
				
				for al in range(2):	
					bubble_allele = bubble.allele0 if al==0 else bubble.allele1
					bubble_allele_name = bubble_allele.name
				
					haplo = "none"
				
					if pred_haplo=="0":
						haplo = "H1" if al==0 else "H2"
					elif pred_haplo == "1":
						haplo = "H2" if al==0 else "H1"
					
					print(bubble_allele_name + '\t' + haplo, file=out)
					