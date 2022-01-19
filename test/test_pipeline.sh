#!/usr/bin/env bash

# run all tests, exit with non-zero code on the first failure
# run from within the test folder

set -eu

trap 'echo "Errorcode $? on line $LINENO" ; exit 1' ERR

export KMERZDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && cd .. && pwd )
export PATH=${KMERZDIR}/scripts:${PATH}

BASENAME=10k

#create a simple random "genome" sequence
random_sequence.py seq 1 10000 0 > ${BASENAME}.flat

#simulate reads from the genome
generate_reads.py ${BASENAME}.flat 20000 100 0.01 0.01 0.5 > ${BASENAME}_reads.flat

#count kmers in the reads
count_kmers.py ${BASENAME}_reads.flat 23 > ${BASENAME}_counts

#assemble the kmers into contigs
assemble_kmers.py ${BASENAME}_counts 30 > ${BASENAME}_assembly.flat

#generate a dot plot comparing kmer positions in original and assembled sequence
compare_kmers.py ${BASENAME}_assembly.flat ${BASENAME}.flat 10 > ${BASENAME}_dots
