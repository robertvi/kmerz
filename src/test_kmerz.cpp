#include "kmerz.h"

#include <cassert>
#include <iostream>

using namespace kmerz;

void test_generate_paths()
{
    //std::string test_genome = "GTCAACGCTACAACTCAAGGGTAAAGCCGAATACACAGACACGGTGTGG"
    //std::string test_genome = "GTCAACGCTACAACTCAAGGGTAAAGCCGAATACACAGACACGGTGTGG"

    std::string kmer1="GTCAACGCTACAACTCAAGGGTAAAGCCGAA";
    std::string kmer2="TCAACGCTACAACTCAAGGGTAAAGCCGAAT";
    std::string kmer3="CAACGCTACAACTCAAGGGTAAAGCCGAATA";
                     //GAACGCTACAACTCAAGGGTAAAGCCGAATA
    std::vector< std::string > test_kmers;

    test_kmers.push_back(kmer1);
    test_kmers.push_back(kmer2);
    test_kmers.push_back(kmer3);

    kmerz::EulerGraph graph(test_kmers);

    graph.generateGraph();

    std::vector< std::string > contig_list;

    graph.generatePaths(contig_list);

    for(auto it=contig_list.begin(); it!=contig_list.end(); it++)
    {
        std::cout << *it << std::endl;
    }
}

void test_reverse_complement()
{
    //assert(reverseComplement("A") == "A");

    assert(reverseComplement("A") == "T");
    assert(reverseComplement("AA") == "TT");
    assert(reverseComplement("AAA") == "TTT");
    assert(reverseComplement("TTTT") == "AAAA");
    assert(reverseComplement("CC") == "GG");
    assert(reverseComplement("G") == "C");

    assert(reverseComplement("AAGGTTCC") == "GGAACCTT");
    assert(reverseComplement("TGCATGCA") == "TGCATGCA");

}

int main(int argc,char*argv[])
{
    test_reverse_complement();
    test_generate_paths();
    return 0;
}
