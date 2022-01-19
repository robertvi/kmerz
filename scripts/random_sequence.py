#!/usr/bin/env python3

'''
generate random sequences in two column flat format
sequential identifiers
controllable proportion of Ns
'''

import sys
import random

#provide usage information
if sys.argv[1] in ['-h','-H','--help','-?']:
    print("Usage: random_sequence.py <identifier> <no_of_sequences> <no_of_bases> [<prob_of_N>]")
    print("results are printed to stdout")
    exit()

identifier = sys.argv[1]
sequences = int(sys.argv[2])
bases = int(sys.argv[3])
prob_N = 0.0

#override default of having no Ns
if len(sys.argv) == 5: prob_N = float(sys.argv[4])

#ensure results are different each time program is run
random.seed()

for id_count in range(sequences):
    #generate bases
    seq = []
    for i in range(bases):
        if random.random() < prob_N:
            seq.append('N')
        else:
            seq.append(random.choice('ATCG'))

    print(identifier + str(id_count) + ' ' + ''.join(seq))
