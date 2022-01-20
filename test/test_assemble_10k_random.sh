#!/usr/bin/env bash
set -eu

BASENAME=10k

mkdir -p tmp/random_10k && cd tmp/random_10k && pwd

#create a simple random "genome" sequence
random_sequence.py seq 1 10000 0 > ${BASENAME}.flat

#simulate reads from the genome
generate_reads.py ${BASENAME}.flat 20000 100 0.01 0.01 0.5 > ${BASENAME}_reads.flat

#count kmers in the reads
count_kmers.py ${BASENAME}_reads.flat 23 > ${BASENAME}_counts

#assemble the kmers into contigs
assemble_kmers.py --kmer-counts=${BASENAME}_counts --min-count=30 > ${BASENAME}_assembly.flat

#generate a dot plot comparing kmer positions in original and assembled sequence
compare_kmers.py ${BASENAME}_assembly.flat ${BASENAME}.flat 10 > ${BASENAME}_dots

cd - > /dev/null
