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
            EulerEdge(uint64_t,EulerNode*,EulerNode*,bool,bool);
            EulerNode*getPrefixNode(){return prefix_node;};
            EulerNode*getSuffixNode(){return suffix_node;};
            bool prefixNodeRevcmp(){return prefix_revcmp;};
            bool suffixNodeRevcmp(){return suffix_revcmp;};
            void setVisited(){visited = true;};
            bool getVisited(){return visited;};

            //get sequence preadjusted for reverse complement status
            //std::string getPrefixSequence();
            std::string getMiddleSequence();
            uint64_t getSequence(){return sequence;};
            char getFirstPrefixBase();
            char getLastSuffixBase();
        private:
            EulerNode*prefix_node; //node containing the first K-1 bases of this kmer
            EulerNode*suffix_node; //node containing the last K-1 bases of this kmer

            //true if node's sequence must be reverse complemented to recover
            //original kmer sequence in it's canonical form
            bool prefix_revcmp;
            bool suffix_revcmp;

            bool visited; //true if this edge was visited already during graph walk

            uint64_t sequence; //canonical form of edge's sequence ie the kmer
    };

    class EulerNode
    {
        public:
            //create node from canonical prefix/suffix sequence
            EulerNode(uint64_t psfix);
            void addIncoming(EulerEdge*);
            void addOutgoing(EulerEdge*);
            uint64_t getSequence(){return sequence;};
            EulerEdge*unvisitedOutgoing();
            EulerEdge*unvisitedIncoming();
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

            //load kmer sequences from in memory vector of strings
            EulerGraph(std::vector< std::string >&);

            //load all kmers into the graph sequentially
            void generateGraph();

            //trace path(s) through all edges until no more left
            void generatePaths(std::vector< std::string > &contig_list);
        private:
            //lookup or generate a node
            void generateNode(uint64_t psfix,EulerNode*&node,bool&revcmp);
            std::string walkPath(EulerEdge*,bool);

            //canonical kmers
            std::vector< uint64_t > kmer_list;

            //graph edges
            std::list< EulerEdge* > graph_edges;

            //graph nodes keyed by their canonical k-1 base sequence
            std::unordered_map< uint64_t, EulerNode* > graph_nodes;
    };

    //convert any length string to uint64_t
    //do not make canonical
    uint64_t anyString2uint64t(const std::string&seq);

    //convert string of bases to bits in a uint64_t
    uint64_t string2uint64t(const std::string&);

    //convert uint64_t kmer into string
    std::string uint64t2string(uint64_t kmer);

    //convert the prefix/suffix k-1 mer into canonical form
    //if not already, return true if change was made
    bool psfix2canonical(uint64_t*psfix);

    //return reverse complement of the k-1 mer
    uint64_t revcmpPsfix(uint64_t psfix);

    //convert k-1 mer to string
    std::string psfix2string(uint64_t psfix);

    //reverse complement a string
    std::string reverseComplement(const std::string&);

    uint64_t kmer2prefix(uint64_t);
    uint64_t kmer2suffix(uint64_t);
} //namespace kmerz
#endif //__ROBERTVI_KMERZ_KMERZ_H__
