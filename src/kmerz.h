#ifndef __ROBERTVI_KMERZ_KMERZ_H__
#define __ROBERTVI_KMERZ_KMERZ_H__

#include <string>
#include <vector>

//#define KMERSIZE 31

namespace kmerz
{
    class EulerGraph
    {
        public:
            //load kmer sequences and counts from file
            EulerGraph(const std::string&inputFile,int minKmerCount);
        private:
            //simple list of canonical kmers
            std::vector< uint64_t > kmer_list;
    };

    uint64_t string2uint64t(const std::string&);
} //namespace kmerz
#endif //__ROBERTVI_KMERZ_KMERZ_H__
