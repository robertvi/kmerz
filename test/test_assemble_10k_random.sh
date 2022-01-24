#!/usr/bin/env bash
set -eu

BASENAME=10k

mkdir -p tmp/random_10k && cd tmp/random_10k && pwd

#create a simple random "genome" sequence of 10000 bases
random_sequence.py seq 1 10000 0 > ${BASENAME}.flat

#simulate reads from the genome
#20000 100 base reads with 1% read errors 1% missing values and 50% from the reverse strand
generate_reads.py ${BASENAME}.flat 20000 100 0.01 0.01 0.5 > ${BASENAME}_reads.flat

#count 23mers in the reads
count_kmers.py ${BASENAME}_reads.flat 23 > ${BASENAME}_counts

#assemble the kmers into contigs filtering out kmers with a count below 30
assemble_kmers.py --kmer-counts=${BASENAME}_counts --min-count=30 > ${BASENAME}_assembly.flat

#generate a dot plot comparing kmer positions in original and assembled sequence
#compare using 10mer
compare_kmers.py ${BASENAME}_assembly.flat ${BASENAME}.flat 10 > ${BASENAME}_dots

cd - > /dev/null
