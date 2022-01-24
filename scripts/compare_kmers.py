#!/usr/bin/env python3

'''
compare kmer positions between two sets of sequences
output coordinates of matching kmers suitable for drawing a dot plot
'''

import sys
from kmer_module import *
from collections import defaultdict

#provide usage information
if sys.argv[1] in ['-h','-H','--help','-?']:
    print("Usage: compare_kmers.py <flat_file1> <flat_file2> <kmer_size>")
    print("results are printed to stdout")
    exit()

input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
kmersize = int(sys.argv[3])

kmer_posn = defaultdict(list)

#index position of all kmers in first sequence set
posn = 0
for identifier,sequence,quality in generate_sequences(input_file1):
	for i,kmer in enumerate(generate_kmers(sequence,kmersize)):
		rev = reverse_complement(kmer)
		if rev < kmer: kmer = rev
		kmer_posn[kmer].append(posn+i)

	posn += len(sequence)

#scan second sequence set looking up in the first set
posn = 0
for identifier,sequence,quality in generate_sequences(input_file2):
	for i,kmer in enumerate(generate_kmers(sequence,kmersize)):
		rev = reverse_complement(kmer)
		if rev < kmer: kmer = rev
		if not kmer in kmer_posn: continue

		for x in kmer_posn[kmer]:
			print(str(x)+' '+str(posn+i))

	posn += len(sequence)
