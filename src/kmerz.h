#ifndef __ROBERTVI_KMERZ_KMERZ_H__
#define __ROBERTVI_KMERZ_KMERZ_H__

#include <string>
#include <vector>
#include <list>

namespace kmerz
{
    class EulerNode;
    
    class EulerEdge
    {
        public:

        private:
            EulerNode*prefix;
            EulerNode*suffix;

            //true if node's sequence must be reverse complemented to recover
            //original kmer sequence as it's canonical form
            bool prefix_reversed;
            bool suffix_reversed;
            bool visited; //true if visited already during graph walk
    };

    class EulerNode
    {
        public:

        private:
            std::list<EulerEdge*>incoming; //edges with this node as their suffix
            std::list<EulerEdge*>outgoing; //edges with this node as their prefix
            uint64_t sequence; //canonical form of node's sequence
    };

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
