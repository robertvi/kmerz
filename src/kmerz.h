#ifndef __ROBERTVI_KMERZ_KMERZ_H__
#define __ROBERTVI_KMERZ_KMERZ_H__

#include <string>
#include <vector>
#include <list>
#include <unordered_map>

namespace kmerz
{
    class EulerNode;

    class EulerEdge
    {
        public:

        private:
            EulerNode*prefix; //node containing the first K-1 bases of this kmer
            EulerNode*suffix; //node containing the last K-1 bases of this kmer

            //true if node's sequence must be reverse complemented to recover
            //original kmer sequence in it's canonical form
            bool prefix_reversed;
            bool suffix_reversed;

            bool visited; //true if this edge was visited already during graph walk
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
            //load kmer sequences into kmer_list with min-count filtering
            EulerGraph(const std::string&inputFile,int minKmerCount);

            //load all kmers into the graph sequentially
            void generateGraph();
        private:
            //simple list of canonical kmers before they go into the graph
            std::vector< uint64_t > kmer_list;

            //graph nodes keyed by their canonical sequence
            std::unordered_map< uint64_t, EulerNode* > graph_nodes;
    };

    uint64_t string2uint64t(const std::string&);
} //namespace kmerz
#endif //__ROBERTVI_KMERZ_KMERZ_H__
