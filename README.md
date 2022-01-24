# kmerz
kmerz is a very simple de Bruijn graph denovo assembly program written in C++ with associated tools in Python. kmerz is not a practical tool to assemble a genome but rather a simple demonstration program to illustrate how de Bruijn graph assembly works and to provide a context for developing other simple demonstration tools for training purposes.

# What is kmerz?
kmerz is currently able to assemble simulated reads from simple randomly generated DNA sequences back into contigs. The scripts folder contains Python tools to generate random "genome" sequences, simulated reads and for counting kmers in the reads.

The intension is to gradually extend the program to include more realistic operations dealing with read errors, heterozygous genomes and repetitive elements etc

## Requirements

C++11 compiler (tested with g++ 9.3.0) and Python 3.
Optional: make, R, ggplot2

## Installation

    #clone the repository
    git clone https://github.com/robertvi/kmerz.git
    cd kmerz

    #build the main binary
    make

    #or to also build and run the tests of the main binary
    ./build_and_test.sh

## Testing

    #to test the Python utilities
    ./test/test_python_scripts.sh

    #to assemble a simple test genome using kmerz
    ./test/test_kmers.sh
