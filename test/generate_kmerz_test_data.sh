#!/usr/bin/env bash
set -eu

BASENAME=10k

mkdir -p tmp/kmerz_test_data && cd tmp/kmerz_test_data && pwd

#create a simple random "genome" sequence
random_sequence.py seq 1 10000 0 > ${BASENAME}.flat

#simulate reads from the genome, without errors of any kind
generate_reads.py ${BASENAME}.flat 20000 100 0.0 0.0 0.5 > ${BASENAME}_reads.flat

#count kmers in the reads
count_kmers.py ${BASENAME}_reads.flat 31 > ${BASENAME}_31mer_counts

#plot histogram to Rplots.pdf
plot_kmer_histogram.R

cd - > /dev/null
