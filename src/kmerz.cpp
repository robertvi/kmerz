#include "kmerz.h"

#include <iostream>
#include <fstream>

namespace kmerz
{

EulerGraph::EulerGraph(const std::string&inputFile,int minKmerCount)
{
    std::ifstream ifs;
    ifs.open (inputFile, std::ios::binary);

    std::string fwd_kmer,rev_kmer;
    int line_ctr=0,count;

    while(ifs.good())
    {
        line_ctr++;
        ifs >> fwd_kmer >> rev_kmer >> count;

        if(fwd_kmer.length() != 31)
        {
            throw std::runtime_error("Error: kmer not 31 bases on line: " + std::to_string(line_ctr) + ":"+fwd_kmer);
        }

        std::cout << fwd_kmer << " " << rev_kmer << " " << count << std::endl;
    }

    ifs.close();
}


} //namespace kmerz
