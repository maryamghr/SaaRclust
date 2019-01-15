import fileinput

bubblenum = 0
allelenum = 1
bubblelines = ''
chainlens = []
allelelentoline = {}
prevline = ''

def printalleles(allelelentoline):
	for allelelen in allelelentoline:
		if len(allelelentoline[allelelen]) > 1:
			for allele in allelelentoline[allelelen]:
				print(allele)

for l in fileinput.input():
	if l[0] == '>':
		prevline = l
		name = l.split('_')
		#allelenum = int(name[3])
		if int(name[1]) != bubblenum:
			# new bubble, update the bubble number and process the previous bubble
			bubblenum = int(name[1])
			printalleles(allelelentoline)
			allelelentoline = {}

	else:
		# still in the same bubble, update allelelentoline
		if len(l) in allelelentoline:
			allelelentoline[len(l)].add(prevline + l.strip())
		else:
			allelelentoline[len(l)] = set([prevline + l.strip()])

printalleles(allelelentoline)
