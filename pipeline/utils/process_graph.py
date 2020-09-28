import pysam
import time

def get_node_haplotypes(haplo_file):
	node_to_haplo = {}
	with open(haplo_file) as f:
		# skip header
		next(f)
		for line in f:
			sp = line.split()
			node_to_haplo[sp[0]] = sp[-1]
	
	return node_to_haplo

def separate_haplotypes(node_to_haplo, input_gfa, h1_gfa, h2_gfa):
	read_haplo = {}
	
	with open(input_gfa) as gfa:
		with open(h1_gfa, 'w') as out1:
			with open(h2_gfa, 'w') as out2:
				
				for line in gfa:
					only_h1, only_h2 = False, False
					sp = line.split()
					node = sp[1]
					
					if not node in node_to_haplo:
						# unitig does not have any ss unique match
						node_to_haplo[node] = "none"
						
					if sp[0] == 'A':
						# the line containing ccs read info
						read_haplo[sp[4]] = node_to_haplo[node]
											
					if sp[0]!="L":
						# non-edge line
						if node_to_haplo[node]=="H1":
							only_h1=True
						elif node_to_haplo[node]=="H2":
							only_h2=True							
						
					else:
						# edge-line
						out_node = sp[3]
						
						edge_haplo = [node_to_haplo[node], node_to_haplo[out_node]]
						
						if "H1" in edge_haplo and "H2" in edge_haplo:
							# edge with endpoints from two different haplotypes; skip
							continue
							
						if "H1" in edge_haplo:
							only_h1=True
						elif "H2" in edge_haplo:
							only_h2=True
							
					
					if not only_h1 and not only_h2:
						print(line.strip(), file=out1)
						print(line.strip(), file=out2)
						continue
						
					if only_h1:
						print(line.strip(), file=out1)
						continue
						
					if only_h2:
						print(line.strip(), file=out2)
						
	return read_haplo
	
def output_read_haplo(read_haplo, fasta_file, output_file):
	'''
	'''

	start_time = time.time()
	print('getting long reads')

	long_reads = {}	
	
	print('getting long reads from', fasta_file, '...')
	
	with open(output_file, 'w') as out:
		with pysam.FastxFile(fasta_file) as fastafile:
			for read in fastafile:
				read_name = read.name
			
				if read_name not in read_haplo:
					print(read_name + "\tnone", file=out)
					continue
					
				print(read_name + "\t" + read_haplo[read_name], file=out)
				
	print('elapsed time =', time.time()-start_time)
	