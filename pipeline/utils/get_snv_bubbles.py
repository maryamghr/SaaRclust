import fileinput

k = int(snakemake.wildcards['k'])
bubblenum = 0
allelenum = 1
bubblelines = ''
chainlens = []
allelelentoline = {}
prevline = ''

revcomp = {'a':'t', 'c':'g', 'g':'c', 't':'a', 'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def reversecomp(seq):
	rc = ''
	for i in range(len(seq)):
		rc = revcomp[seq[i]] + rc
	return rc


def hamming_dist(a, b):
	if (len(a) != len(b)):
		print("the two sequences are not of the same length")
		return None
	return sum([1 if a[i]!=b[i] else 0 for i in range(len(a))])


def printalleles(allelelentoline, k, out):
	for allelelen in allelelentoline:
		# consider only biallelic bubbles with chains of the same length
		if len(allelelentoline[allelelen]) == 2:
			# consider only bubbles with reasonably low hamming distance in their chains as valid snv bubbles
			valid = False
			s0 = allelelentoline[allelelen][0].split("\n")[-1]
			s1 = allelelentoline[allelelen][1].split("\n")[-1]
			if hamming_dist(s0, s1)*k < len(s0):
				valid = True
			else:
				s1_rc = reversecomp(s1)
				if hamming_dist(s0, s1_rc)*k < len(s0):
					allelelentoline[allelelen][1] = allelelentoline[allelelen][1].split("\n")[0] + "\n" + s1_rc
					valid = True
			
			if valid:
				for allele in allelelentoline[allelelen]:
					print(allele, file=out)


with open(snakemake.input[0]) as f:
	with open(snakemake.output[0], 'w') as out:
		for l in f:
			if l[0] == '>':
				prevline = l
				name = l.split('_')
		
				if int(name[1]) != bubblenum:
					# new bubble, update the bubble number and process the previous bubble
					bubblenum = int(name[1])
					if len(allelelentoline) > 0:
						printalleles(allelelentoline, k, out)
						allelelentoline = {}

			else:
				# update allelelentoline
				if len(l) in allelelentoline:
					allelelentoline[len(l)].append(prevline + l.strip())
				else:
					allelelentoline[len(l)] = [prevline + l.strip()]

		printalleles(allelelentoline, k, out)
