#!/usr/bin/env python3

'''
kmer counting tool
counts canonical kmers in flat format file
simple example tool for illustrative purposes only
'''

import sys
from kmer_module import *
from collections import defaultdict

#provide usage information
if sys.argv[1] in ['-h','-H','--help','-?']:
    print("Usage: count_kmers.py <flat_input_filename> <kmer_size>")
    exit()

input_file = sys.argv[1]
kmer_size = int(sys.argv[2])
kmer_count = defaultdict(int)

with open(input_file) as f:
    for line in f:
        column = line.strip().split()
        identifier = column[0]
        sequence = column[1]

        for kmer in generate_kmers(sequence,kmer_size):
            if 'N' in kmer: continue
            rev = reverse_complement(kmer)
            if rev < kmer: kmer = rev
            kmer_count[kmer] += 1

for kmer in kmer_count:
    print(kmer + ' ' + str(kmer_count[kmer]))
