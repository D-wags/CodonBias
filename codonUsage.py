import sys
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from decimal import *


# reads in file from command line
def readFile(file):
	inFile = open(file, 'r')
	contents = inFile.read()
	inFile.close()
	return contents


# generates dictionary containing every possible codon
def generate_codons():
	bases = ['a', 'c', 'g', 't']
	every_codon = {}
	
	codons = [nuc1+nuc2+nuc3 for nuc1 in bases for nuc2 in bases for nuc3 in bases]
	for codon in codons:
		every_codon[codon] = 0
	return every_codon


# returns list of codons and codon freqs in dna argument
def codonUsage(dna):

	codon_freq = {}
	dna = dna.lower()

	for nuc in range(len(dna)):
		if dna[nuc: nuc+3] == "atg":
			start = dna[nuc: nuc+3]
			orf = dna[nuc:]
			break

	orf = orf.replace('\n', '').replace('\r', '').replace(' ', '')
	length = len(orf) -2

	total = 0
	codon_list = []

	for x in range(0, length, 3):
		total += 1
		codon = orf[x:x+3]

		
		codon_list.append(codon)
		if codon in codon_freq:
			codon_freq[codon] += 1
		else:
			codon_freq[codon] = 1
	print(len(codon_list))


	all_codons = generate_codons()
	for codon in all_codons:
		for gene_codon in codon_freq:
			if codon == gene_codon:
				all_codons[codon] += codon_freq[gene_codon]
	
	codons = []
	codon_freqs = []
	indexes = []

	n = 1
	for codon in all_codons:
		codons.append(codon)
		codon_freqs.append((all_codons[codon]/total)*100)
		indexes.append(n)
		n += 1


	return codons, codon_freqs, indexes, codon_list


# translates codon list and returns list of 20 aa's and list of frequency indexed together
def translate(codons):
	genetic_code = {'phe': ['ttt', 'ttc'], 'leu': ['tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'], 'ile': ['att', 'atc', 'ata'], 
	'met': ['atg'], 'val': ['gtt', 'gtc', 'gta', 'gtg'], 'ser': ['tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'], 'pro': ['cct', 'ccc', 'cca', 'ccg'],
	'thr': ['acc', 'act', 'acg', 'aca'], 'ala': ['gct', 'gcc', 'gca', 'gcg'], 'tyr': ['tat', 'tac'], 'stop': ['taa', 'tag', 'tga'], 
	'his': ['cat', 'cac'], 'gln': ['caa', 'cag'], 'asn': ['aat', 'aac'], 'lys': ['aaa', 'aag'], 'asp': ['gat', 'gac'], 'glu': ['gaa', 'gag'],
	'cys': ['tgt', 'tgc'], 'trp': ['tgg'], 'arg': ['cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'], 'gly': ['ggt', 'gga', 'ggc', 'ggg']  } 

	amino_acids = {"ala": 0, 'pro': 0, "leu": 0, "ile": 0, "trp": 0, 'tyr':0, 'lys': 0, "ser": 0, 'gly': 0, 'cys': 0, 'asn': 0, 'gln': 0, 'asp': 0,
	'glu': 0, 'his':0, 'phe': 0, 'met': 0, 'val': 0, 'arg': 0, "stop": 0, "thr": 0}


	for codon in codons:
		for aa in genetic_code:
			if codon in genetic_code[aa]:
				amino_acids[aa] += 1

	aas = []
	aa_freq = []
	totalpcnt = 0

	for aa in amino_acids:
		aas.append(aa)
		aa_freq.append(Decimal(amino_acids[aa]/len(codons)))
		totalpcnt += amino_acids[aa]/len(codons)
	# print(totalpcnt)
	# print(len(amino_acids))
	#print("Total percent = %d" % totalpcnt)	
	return aas, aa_freq



# plots codon bar chart and amino acid composition pie chart
def plotter(size, cfreqs, codons, cfreqs2, aa_freq, aas):
	N = len(size)

	ind = np.arange(N)
	width = 0.55
	fig = plt.figure()
	ax = fig.add_subplot(211)
	target = ax.bar(ind, cfreqs, width, color='b',align='center')
	e_coli = ax.bar(ind+width, cfreqs2, width, color='y',align='center')
	

	ax.set_ylabel("Frequency (%)")
	ax.set_xlabel("Codons")
	ax.set_title("Codon Usage in Target Gene vs E.coli")
	ax.set_xticks(ind + width / 2)

	ax.set_xticklabels(codons, rotation='vertical')
	ax.legend((target[0], e_coli[0]), ('Target', 'E.coli'))

	ax1 = fig.add_subplot(212)
	ax1.pie(aa_freq, labels=aas, autopct='%1.1f%%', shadow=False, startangle=90)
	ax1.axis('equal')

	plt.show()

gene1 = readFile(sys.argv[1])
gene2 = readFile(sys.argv[2])

#target gene
codons, codon_freqs, indexes, codon_list = codonUsage(gene1)

# ecoli reference gene
codons1, codon_freqs1, indexes1, codon_list1 = codonUsage(gene2)

trgt_aas, trgt_freq = translate(codon_list)


plotter(indexes, codon_freqs, codons, codon_freqs1, trgt_freq, trgt_aas)








