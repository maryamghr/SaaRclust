#!/usr/bin/python

import fileinput

edges = {}
nodes = {}

rev_end = {"+": "-", "-": "+"}

def canon_edgepair(left, right):
	if right[0] < left[0]: return (right, left)
	elif left[0] < right[0]: return (left, right)

def canon_edge(left, right):
	if right[0] < left[0]: return ((right[0], rev_end[right[1]]), (left[0], rev_end[left[1]]))
	elif left[0] < right[0]: return (left, right)

for l in fileinput.input():
	parts = l.split('\t')
	if parts[0] == 'L':
		fw_source = (parts[1], parts[2])
		fw_target = (parts[3], parts[4])
		bw_source = (parts[3], rev_end[parts[4]])
		bw_target = (parts[1], rev_end[parts[2]])
		if fw_source not in edges: edges[fw_source] = set()
		if bw_source not in edges: edges[bw_source] = set()
		edges[fw_source].add(fw_target)
		edges[bw_source].add(bw_target)
	if parts[0] == 'S':
		nodes[parts[1]] = parts[2].strip()

alleles = {}

for node in nodes:
	fw = (node, "+")
	bw = (node, "-")
	if fw not in edges: continue
	if bw not in edges: continue
	if len(edges[fw]) != 1: continue
	if len(edges[bw]) != 1: continue
	left = edges[bw].pop()
	right = edges[fw].pop()
	if left[0] == right[0]: continue
	canon = canon_edgepair(left, right)
	if canon not in alleles: alleles[canon] = []
	alleles[canon].append(node)

bubblenum = 1

for bubble in alleles:
	if len(alleles[bubble]) == 1: continue
	#print bubbles in whatever format you want. alleles[bubble] contains all nodes which are a part of this bubble
	allelenum = 1
	for node in alleles[bubble]:
		print(">bubble_" + str(bubblenum) + "_allele_" + str(allelenum) + "_node_" + str(node))
		allelenum += 1
		print(nodes[node])
	bubblenum += 1
