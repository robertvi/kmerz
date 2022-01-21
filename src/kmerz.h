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
            EulerEdge(EulerNode*,EulerNode*,bool,bool);
        private:
            EulerNode*prefix_node; //node containing the first K-1 bases of this kmer
            EulerNode*suffix_node; //node containing the last K-1 bases of this kmer

            //true if node's sequence must be reverse complemented to recover
            //original kmer sequence in it's canonical form
            bool prefix_revcmp;
            bool suffix_revcmp;

            bool visited; //true if this edge was visited already during graph walk
    };

    class EulerNode
    {
        public:
            //create node from canonical prefix/suffix sequence
            EulerNode(uint64_t psfix);
            void addIncoming(EulerEdge*);
            void addOutgoing(EulerEdge*);
        private:
            std::list< EulerEdge* > incoming; //edges with this node as their suffix
            std::list< EulerEdge* > outgoing; //edges with this node as their prefix
            uint64_t sequence; //canonical form of node's sequence
    };

    class EulerGraph
    {
        public:
            //load kmer sequences into kmer_list with min-count filtering
            EulerGraph(const std::string&inputFile,int minKmerCount);

            //load all kmers into the graph sequentially
            void generateGraph();

            //trace path(s) through all edges until no more left
            void generatePaths(const std::string&);
        private:
            //lookup or generate a node
            void generateNode(uint64_t psfix,EulerNode*&node,bool&revcmp);

            //simple list of canonical kmers before they go into the graph
            std::vector< uint64_t > kmer_list;

            //graph edges
            std::list< EulerEdge* > graph_edges;

            //graph nodes keyed by their canonical k-1 base sequence
            std::unordered_map< uint64_t, EulerNode* > graph_nodes;
    };

    //convert string of bases to bits in a uint64_t
    uint64_t string2uint64t(const std::string&);

    //convert the prefix/suffix k-1 mer into canonical form
    //if not already, return true if change was made
    bool psfix2canonical(uint64_t*psfix);
} //namespace kmerz
#endif //__ROBERTVI_KMERZ_KMERZ_H__
