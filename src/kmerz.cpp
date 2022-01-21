#include "kmerz.h"

#include <iostream>
#include <fstream>

namespace kmerz
{
const uint64_t KMER_SIZE=31;
const uint64_t KMER_MASK=( ((uint64_t)1)<<(KMER_SIZE*2) )-1;
const uint64_t KMER_FIRSTOFFSET=(KMER_SIZE-1)*2;

EulerGraph::EulerGraph(const std::string&inputFile,int minKmerCount)
{
    std::ifstream ifs;
    ifs.open (inputFile, std::ios::binary);

    std::string fwd_kmer,rev_kmer;
    int line_ctr=0;
    uint64_t count;

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

        if(count < minKmerCount) continue;

        //convert kmer to uint64_t representation
        uint64_t canonical_kmer = string2uint64t(fwd_kmer.c_str());
        kmer_list.push_back( canonical_kmer );

        std::cout << fwd_kmer << " " << canonical_kmer << " " << count << std::endl;
    }

    ifs.close();
}

//converts to uint64_t and ensures its canonical
uint64_t string2uint64t(const std::string&seq)
{
    uint64_t fkmer=0,rkmer=0;

    for(auto p=0; p<KMER_SIZE; p++)
    {
        //fkmer: grows from the right (LSB)
        //rkmer: grows from the left (MSB)
        switch(seq[p])
        {
            case 'A':
            case 'a':
                fkmer = ((fkmer << 2) +  uint64_t(0)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(3) << KMER_FIRSTOFFSET);
                break;
            case 'C':
            case 'c':
                fkmer = ((fkmer << 2) +  uint64_t(1)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(2) << KMER_FIRSTOFFSET);
                break;
            case 'G':
            case 'g':
                fkmer = ((fkmer << 2) +  uint64_t(2)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(1) << KMER_FIRSTOFFSET);
                break;
            case 'T':
            case 't':
                fkmer = ((fkmer << 2) +  uint64_t(3)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(0) << KMER_FIRSTOFFSET);
                break;
            default:
                throw std::runtime_error("Error: invalid character in kmer: " + seq[p]);
        }
    }

    //return canonical version
    if (fkmer <= rkmer) return fkmer;
    else                return rkmer;
}

} //namespace kmerz
