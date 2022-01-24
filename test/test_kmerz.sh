#!/usr/bin/env bash
set -eu

BASENAME=10k

export KMERZDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && cd .. && pwd )
export PATH=${KMERZDIR}/scripts:${PATH}
export PATH=${KMERZDIR}:${PATH}

mkdir -p tmp/kmerz_test_data && cd tmp/kmerz_test_data && pwd

#create a simple random "genome" sequence
random_sequence.py seq 1 10000 0 > ${BASENAME}.flat

#simulate 100 base reads from the genome, without errors of any kind
generate_reads.py ${BASENAME}.flat 20000 100 0.0 0.0 0.5 > ${BASENAME}_reads.flat

#count 31mers in the reads
count_kmers.py ${BASENAME}_reads.flat 31 > ${BASENAME}_31mer_counts

#assemble using kmerz
kmerz -i ${BASENAME}_31mer_counts > ${BASENAME}_assembly.flat

#generate coordinates for a dot plot comparing kmer positions in original and assembled sequence
#compare using 20mers
compare_kmers.py ${BASENAME}_assembly.flat ${BASENAME}.flat 20 > ${BASENAME}_dots

#generate the dot plot using ggplot
make_dot_plot.R ${BASENAME}_dots ${BASENAME}_dotplot.png
