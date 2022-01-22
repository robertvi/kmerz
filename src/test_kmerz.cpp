#include "kmerz.h"

#include <cassert>
#include <iostream>

using namespace kmerz;

std::vector<std::string> test_list =
{
"CAGCATTATGAATGAACGAAAGTGTGGGCCG",
"TCCTAGTTCTTTGTCTCGACATAATGCGACG",
"AAGAACCGGCTATCTTGAGCTGAAGGCTGGC",
"TCTTTCGCTACAGGAATACCATTAAAAAGCG",
"TCCGTAGATCGAAGATCTCCCCCAAACTCGG",
"TGGGAGAGAGTCATATAGCGTTTGTTGTAAC",
"GCCCGCGTCCTGAGTTGAAAAGGACACTGCA",
"TACGAGGATTAATAGATTATTTCCAAGAAAT",
"GGGTGCTAACTTGCGGATGGTTGAGACTATG",
"TGGGCTAATGTAAGCCCATCATGACATTTCT",
"CCCCTAAAGTACTCAGAACTTGTGTACTCGC",
"TACGTAAGTCAGACGTCGGACGATTAACTCT",
"GCCGATAAACGAAAACTCATGTGTCCACTAA",
"CCCGTAATACGGCGGGCGCCGGGCAAGCCGG",
"CATGGATATATAGGTAAGAAGCTCTCTACTG",
"CTTACATATGCCAGAGCGCGTATCCCGAGGT",
"CTGATAGAATGTGAGTGGAGGGAAAACAGTA",
"GGGATGGTACGTCATAGTGATACGGCGACAT",
"GCCTTCGGGAACGGTGCACACAGCGGGCACC",
"GGCATTCAGCAGCTATTACTAGACATGATAG",
"CACGATAAGAAACGGCCGCCATGGAAAATGC"
};

void test_utility_functions()
{
    //assert(reverseComplement("A") == "XYZ");

    assert(reverseComplement("A") == "T");
    assert(reverseComplement("AA") == "TT");
    assert(reverseComplement("AAA") == "TTT");
    assert(reverseComplement("TTTT") == "AAAA");
    assert(reverseComplement("CTTT") == "AAAG");
    assert(reverseComplement("ATTT") == "AAAT");
    assert(reverseComplement("GTTT") == "AAAC");
    assert(reverseComplement("CC") == "GG");
    assert(reverseComplement("G") == "C");
    assert(reverseComplement("AAGGTTCC") == "GGAACCTT");
    assert(reverseComplement("TGCATGCA") == "TGCATGCA");

    for(int i=0; i<test_list.size(); i++)
    {
        std::string seq = test_list[i];

        assert(reverseComplement(reverseComplement(seq)) == seq);
    }

    for(int i=0; i<test_list.size()-1; i++)
    {
        std::string seq1 = test_list[i];
        std::string seq2 = reverseComplement(seq1);

        if(seq1.compare(seq2) > 0)
        {
            assert(makeCanonical(seq1) == true);
        }
        else
        {
            assert(makeCanonical(seq1) == false);
        }
    }
}

int main(int argc,char*argv[])
{
    test_utility_functions();
    //test_generate_paths();
    return 0;
}
