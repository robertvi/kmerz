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

            //walk forward from a given seed kmer
            std::string walkForwards(std::string);

            //report only unique possible extensions of the seed string
            char extendSuffix(std::string);

            //print as FASTA format
            void printToStdout(std::vector< std::string >&,int=80);

            //print in flat format, one sequence per line
            void printToStdoutFlat(std::vector< std::string >&);
        private:

            //canonical kmers as strings
            std::unordered_set< std::string > kmer_set;
    };

    //reverse complement a string
    std::string reverseComplement(const std::string&);

    //if not canonical reverse complement and return true
    bool makeCanonical(std::string&);
} //namespace kmerz
#endif //__ROBERTVI_KMERZ_KMERZ_H__
