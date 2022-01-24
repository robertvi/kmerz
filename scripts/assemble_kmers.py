#!/usr/bin/env python3

'''
read in kmer counts
assemble into contigs based on simple frequency cut off
'''

import sys
import argparse
from kmer_module import *

def walk_forwards(contig,kmer_length,kmer_counts):
    'walk the contig forwards through all unique extensions'

    while True:
        #extract the last kmer-length-minus-one bases of the contig
        starting_seq = contig[-(kmer_length-1):]

        #search for a unique extension in the remaining kmers
        next_base = extend_sequence(starting_seq,kmer_counts)

        #walk ends when no unique extension is found
        if next_base == None: return contig

        #extend the contig using the unique match
        contig += next_base

def extend_sequence(seq,kmer_counts):
    'try to forward-extend seq by one base using the kmer dictionary'

    next_base = None
    unique = True

    #find all possible extensions and remove them from the kmer dictionary
    for base in 'ATCG':
        kmer = seq + base

        #ensure we only look up canonical kmers
        rev = reverse_complement(kmer)
        if rev < kmer: kmer = rev

        if kmer in kmer_counts:
            del kmer_counts[kmer]
            if next_base == None:
                #possible unique extension
                next_base = base
            else:
                #extension is not unique
                unique = False

    if next_base != None and unique == True: return next_base

    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--min-count', type=int,required=True, help='Minimum kmer count filter')
    parser.add_argument('--kmer-counts', type=str, required=True, help='File containung kmer counts')

    args = parser.parse_args()

    input_file = args.kmer_counts
    min_count = args.min_count

    kmer_counts = {}

    with open(input_file) as f:
        for line in f:
            column = line.strip().split()
            kmer = column[0]
            count = int(column[1])

            if count < min_count: continue

            kmer_counts[kmer] = count

    contig_list = []

    while len(kmer_counts) > 0:
        #pick a seed kmer from the remaining kmers
        seed_kmer = next(iter(kmer_counts))
        kmer_length = len(seed_kmer)
        contig = seed_kmer
        del kmer_counts[seed_kmer]

        #walk the contig forward through the kmer graph
        contig = walk_forwards(contig,kmer_length,kmer_counts)

        #reverse compliment the contig to allow walking backwards
        contig = reverse_complement(contig)

        #walk the (now reversed) contig forward through the kmer graph again
        contig = walk_forwards(contig,kmer_length,kmer_counts)

        #append completed contig to final list
        contig_list.append(contig)

    #print contigs
    for seq_counter,contig in enumerate(contig_list):
        print("contig" + str(seq_counter) + ' ' + contig)
