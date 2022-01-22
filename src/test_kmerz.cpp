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

    /*for(auto it=contig_list.begin(); it!=contig_list.end(); it++)
    {
        std::cout << *it << std::endl;
    }*/
}

void test_utility_functions()
{
    //assert(reverseComplement("A") == "XYZ");

    assert(reverseComplement("A") == "T");
    assert(reverseComplement("AA") == "TT");
    assert(reverseComplement("AAA") == "TTT");
    assert(reverseComplement("TTTT") == "AAAA");
    assert(reverseComplement("CC") == "GG");
    assert(reverseComplement("G") == "C");

    assert(reverseComplement("AAGGTTCC") == "GGAACCTT");
    assert(reverseComplement("TGCATGCA") == "TGCATGCA");

    //assert(anyString2uint64t("AAAAAAAA") == 9999);


    assert(anyString2uint64t("AAAAAAAT") == 3);
    assert(anyString2uint64t("C") == 1);
    assert(anyString2uint64t("G") == 2);
    assert(anyString2uint64t("T") == 3);

    assert(anyString2uint64t("TTGCA") == (3<<8) + (3<<6) + (2<<4) + (1<<2) + (0));

    assert(anyString2uint64t("GCATCG") == (2<<10) + (1<<8) + (0<<6) + (3<<4) + (1<<2) + (2<<0));

    assert( (~(0x1) & 0x3) == 0x2 );

    std::string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    uint64_t intseq = anyString2uint64t(seq);
    assert(intseq == 0);
    std::string revcmpseq = reverseComplement(seq);
    assert(revcmpseq == "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    uint64_t intrevcmpseq = anyString2uint64t(revcmpseq);
    assert(intrevcmpseq == 0xfffffffffffffff);

    assert(revcmpPsfix(intseq) == intrevcmpseq);
    //assert( revcmpPsfix(prefix) == prefix_revcmp );

    for(int i=0; i<test_list.size(); i++)
    {
        //assert( anyString2uint64t(test_list[i]) == string2uint64t(test_list[i])+1 );

        assert( anyString2uint64t(test_list[i]) == string2uint64t(test_list[i]) || \
                anyString2uint64t(reverseComplement(test_list[i])) == string2uint64t(test_list[i]));

        uint64_t kmer = anyString2uint64t(test_list[i]);

        std::cout << test_list[i] << " " << std::hex << anyString2uint64t(test_list[i]) << std::endl;

        std::string prefix_str = test_list[i].substr(0,test_list[i].size()-1);
        std::cout << prefix_str << " " << std::hex << anyString2uint64t(prefix_str) << std::endl;

        std::string prefix_str_revcmp = reverseComplement(prefix_str);
        std::cout << prefix_str_revcmp << std::endl;

        uint64_t prefix = anyString2uint64t(prefix_str);
        uint64_t prefix_copy = prefix;
        uint64_t prefix_revcmp = anyString2uint64t(prefix_str_revcmp);

        assert( kmer2prefix(kmer) == prefix );
        std::cout << std::hex << prefix << " " << revcmpPsfix(prefix) << " " << prefix_revcmp << std::endl;
        //assert( revcmpPsfix(prefix) == prefix_revcmp );

        std::string suffix_str = test_list[i].substr(1,test_list[i].size()-1);
        std::string suffix_str_revcmp = reverseComplement(suffix_str);

        std::cout << test_list[i] << std::endl;
        std::cout << suffix_str << std::endl;
        std::cout << suffix_str_revcmp << std::endl;

        uint64_t suffix = anyString2uint64t(suffix_str);
        uint64_t suffix_copy = suffix;
        uint64_t suffix_revcmp = anyString2uint64t(suffix_str_revcmp);

        std::cout << std::hex << anyString2uint64t(test_list[i]) << std::endl;
        std::cout << std::hex << suffix << std::endl;
        std::cout << std::hex << suffix_revcmp << std::endl;
        std::cout << psfix2string(suffix_revcmp) << std::endl;
        std::cout << suffix_str_revcmp << std::endl;

        //assert( false );
        assert( kmer2suffix(kmer) == suffix );

        std::cout << std::hex << suffix_revcmp << std::endl;
        std::cout << std::hex << revcmpPsfix(suffix) << std::endl;
        assert( revcmpPsfix(suffix) == suffix_revcmp );

        bool result = psfix2canonical(&prefix);
        if(prefix_revcmp < prefix)
        {
            //assert(result == true && prefix == prefix_revcmp);
        }
        else
        {
            //assert(result == false && prefix == prefix_copy);
        }

        result = psfix2canonical(&suffix);
        if(suffix_revcmp < suffix)
        {
            assert(result == true && suffix == suffix_revcmp);
        }
        else
        {
            assert(result == false && suffix == suffix_copy);
        }
    }


}

int main(int argc,char*argv[])
{
    test_utility_functions();
    test_generate_paths();
    return 0;
}
