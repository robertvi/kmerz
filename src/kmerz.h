#ifndef __ROBERTVI_KMERZ_KMERZ_H__
#define __ROBERTVI_KMERZ_KMERZ_H__

#include <string>
#include <vector>
#include <list>
#include <unordered_set>

namespace kmerz
{
    class EulerGraph
    {
        public:
            //load kmer sequences into kmer_list with min-count filtering
            EulerGraph(const std::string&inputFile,int minKmerCount);

            //load kmer sequences from in memory vector of strings
            EulerGraph(std::vector< std::string >&);

            //walk through the graph from seed kmers
            void generateContigs(std::vector< std::string >&);

        private:
            std::string walkForwards(std::string);
            char extendSuffix(std::string);

            //canonical kmers as strings
            std::unordered_set< std::string > kmer_set;
    };

    //return the last size-1 characters
    std::string getSuffix(const std::string&);

    //reverse complement a string
    std::string reverseComplement(const std::string&);

    //if not canonical reverse complement and return true
    bool makeCanonical(std::string&);
} //namespace kmerz
#endif //__ROBERTVI_KMERZ_KMERZ_H__
