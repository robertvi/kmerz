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

        ifs >> fwd_kmer;
        
        if(ifs.eof()) break;

        if(ifs.fail())
        {
            throw std::runtime_error("Error: failed to read forward kmer on line: " + std::to_string(line_ctr));
        }
        if(fwd_kmer.length() != 31)
        {
            throw std::runtime_error("Error: forward kmer not 31 bases on line: " + std::to_string(line_ctr) + ":"+fwd_kmer);
        }

        ifs >> rev_kmer;
        if(ifs.fail())
        {
            throw std::runtime_error("Error: failed to read reverse kmer on line: " + std::to_string(line_ctr));
        }
        if(rev_kmer.length() != 31)
        {
            throw std::runtime_error("Error: reverse kmer not 31 bases on line: " + std::to_string(line_ctr) + ":"+fwd_kmer);
        }

        ifs >> count;
        if(ifs.fail())
        {
            throw std::runtime_error("Error: failed to read kmer count on line: " + std::to_string(line_ctr));
        }
        if(count < 1)
        {
            throw std::runtime_error("Error: invalid kmer count on line: " + std::to_string(line_ctr) + ": count was " + std::to_string(count));
        }

        std::cout << fwd_kmer << " " << rev_kmer << " " << count << std::endl;
    }

    ifs.close();
}


} //namespace kmerz
