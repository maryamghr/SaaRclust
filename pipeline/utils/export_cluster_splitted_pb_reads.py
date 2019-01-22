import sys

clsutpartners = open(sys.argv[0])
pbreads = open(sys.argv[1])


revcomp = {'a':'t', 'c':'g', 'g':'c', 't':'a', 'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def reversecomp(seq):
	rc = ''
	for i in range(len(seq)):
		rc = revcomp[seq[i]] + rc
	return rc


