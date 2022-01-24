# kmerz
kmerz is a very simple de Bruijn graph denovo assembly program written in C++ with associated tools in Python. kmerz is not a practical tool to assemble a genome but rather a simple demonstration program to illustrate how de Bruijn graph assembly works and to provide a context for developing other simple demonstration tools for training purposes.

# What can kmerz do?
kmerz is currently able to assemble simulated reads from simple randomly generated DNA sequences back into contigs. The scripts folder contains Python tools to generate random "genome" sequences, create simulated reads from them and for counting kmers in those reads. The read counts are then fed into kmerz which tries to assemble them.

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

## De Bruijn graph de novo assembly
So called next generation sequencing machines cannot read long stretches of DNA, instead the genome is broken into short fragments which then produce "reads" of 100 to 250 bases. The problem of reassembling the original genome sequence from these short reads is known as de novo genome assembly. The obvious approach of looking for all the possible overlaps between the reads is not feasible due to the large number of reads generated. All the pairwise comparisons would take too long to perform, although this type of approach, called overlap layout consensus genome assembly, was feasible with the first generation of sequencing technologies (Sanger sequencing of DNA cloned into bacteria).

Instead De Bruijn graph assemblers work with kmers. A kmer is a DNA sequence of length k. For example if we pick k=31 then we are considering 31-mers, the default length used by kmerz. De Bruijn graph assemblers work by counting how many times they see each specific kmer sequence in the reads. They then try to reconstruct the original sequence starting from a seed kmer, and seeing which if any of the four possible kmers corresponding to a one base extension (A,T,C or G) of the seed is present in the reads. kmerz, being a very simple assembler, will simply keep extending the sequence as long as it only finds a single possible extension base and will stop if no extensions are found or if more than one possible extension is present, likely indicating the presence of repetitive sequences of some kind. The simple genomes used in the test scripts are purely random and therefore, depending on the kmer and genome lengths, are unlikely to contain any kmers more than once. But real genome contains repetitive sequences of many kinds, and any practical assembler needs a way to deal with them.

Sequencing machines make occassional errors, the most common of which is to misread a single base. If the read depth is great enough then read errors produce kmers which are much rarer than the kmers coming from the genuine genome sequence and can be filtered out using a simple minimum frequency threshold, which kmerz implements with its minimum count parameter.
