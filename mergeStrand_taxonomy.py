#!/usr/bin/env python
import sys

def getValidLineage(nodes):
	lineage = []
	for node in nodes:
		if node.find('unidentified') != -1:
			break
		else:
			lineage.append(node)
	return lineage

def getDict(F):
	Dict = {}
	for line in F[1:]:
		cols = line.strip().split('\t')
		ID = cols[0].strip()
		tax = getValidLineage(cols[1].split(';'))
		conf = float(cols[2])
		Dict[ID] = [tax, conf]
	return Dict

assign1 = getDict(open(sys.argv[1], 'r').readlines())
assign2 = getDict(open(sys.argv[2], 'r').readlines())

print ('Feature ID\tTaxon\tConfidence')
for ID in assign1.keys():
	tax1, conf1 = assign1[ID]
	tax2, conf2 = assign2[ID]
	if (len(tax1) > len(tax2)): 
		out = '\t'.join([ID, ';'.join(tax1), str(conf1)])
	elif (len(tax1) == len(tax2) and conf1 > conf2):
		out = '\t'.join([ID, ';'.join(tax1), str(conf1)])
	else:
		out = '\t'.join([ID, ';'.join(tax2), str(conf2)])
	print (out)
