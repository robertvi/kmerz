#!/usr/bin/env python3

'''
generate simulated reads from sequences in flat format
'''

import sys
import random
from kmer_module import *

#provide usage information
if sys.argv[1] in ['-h','-H','--help','-?']:
    print("Usage: generate_reads.py <input_sequence_filename> <number_of_reads> <read_length> <prob_error> <prob_N> <prob_revcomp>")
    print("input must be in flat format")
    print("results are printed to stdout")
    exit()

input_file = sys.argv[1]
reads = int(sys.argv[2])
length = int(sys.argv[3])
prob_error = float(sys.argv[4])
prob_N = float(sys.argv[5])
prob_revcomp = float(sys.argv[6])

length_check_passed = False

#load in the source sequence(s)
sequence_list = []
with open(input_file) as f:
    for line in f:
        #split line into its columns
        column = line.strip().split()

        #for clarity store each column in a separate variable
        identifier = column[0]
        sequence = column[1]

        #we need at least one sequence longer than the read length
        if len(sequence) >= length:
            length_check_passed = True

        #append to list of sequences
        sequence_list.append([identifier,sequence])

assert(length_check_passed)

#generate the reads
for read_counter in range(reads):
    #pick a random sequence
    sequence_number = random.choice(range(len(sequence_list)))

    #pick a random starting position
    max_position = len(sequence_list[sequence_number][1]) - length
    start = random.randint(0,max_position)

    #extract the read
    read = sequence_list[sequence_number][1][start:start+length]

    #switch 50% to the reverse strand
    if random.random() < prob_revcomp:
        read = reverse_compliment(read)

    #split read into a list of bases for ease of modification
    base_list = [x for x in read]

    #apply errors and Ns
    for i,x in enumerate(base_list):
        if x == 'N': continue

        #convert to N
        if random.random() < prob_N:
            base_list[i] = 'N'

        #introduce a simple read error
        elif random.random() < prob_error:
            j = ('ATCG'.index(x) + random.randint(1,3)) % 4
            base_list[i] = 'ATCG'[j]

    ##convert back to a single string again
    read = ''.join(base_list)

    #print in flat format
    print(sequence_list[sequence_number][0]+'_'+str(read_counter)+' '+read)
