#include "kmerz.h"

#include <iostream>
#include <fstream>

namespace kmerz
{

//how many bases in a kmer, ie what is the value of 'k'
//32 is the maximum possible using a uint64_t
//but k is usually chosen to be odd to prevent the reverse complement of
//the kmer equaling itself
const uint64_t KMER_SIZE=31;

//bit mask set to 1 for all bits used to represent a kmer
const uint64_t KMER_MASK=( ((uint64_t)1)<<(KMER_SIZE*2) )-1;

//bit shift required to shift the two least significant bits
//to the position of the two most significant bits of the kmer
const uint64_t KMER_SHIFT=(KMER_SIZE-1)*2;

//how many bases in the prefix and suffix sequences
const uint64_t PSFIX_SIZE=(KMER_SIZE-1);

//bit mask set to 1 for all bits used to represent a prefix or suffix
//ie a (k-1)mer, *not* the position of a prefix *within* a kmer
const uint64_t PSFIX_MASK=( ((uint64_t)1)<<(PSFIX_SIZE*2) )-1;

//bit shift required to shift the two least significant bits
//to the position of the two most significant bits of the (k-1)mer
const uint64_t PSFIX_SHIFT=(PSFIX_SIZE-1)*2;

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

        //std::cout << fwd_kmer << " " << canonical_kmer << " " << count << std::endl;
    }

    ifs.close();
}

//lookup existing or generate new node for a given prefix or suffix
void EulerGraph::generateNode(uint64_t psfix,bool&revcmp,EulerNode*&node)
{
    //reverse complement the sequence if not already canonical
    revcmp = psfix2canonical(&psfix);

    //lookup existing node with this sequence
    auto node_it = graph_nodes.find(psfix);

    if(node_it == graph_nodes.end())
    {
        //create the missing node
        node = new EulerNode(psfix);

        //store in map
        graph_nodes[psfix] = node;
    }
    else
    {
        //get existing node
        node = node_it->second;
    }
}
//load all kmers into the graph sequentially
//generate the associated mode and edge objects
void EulerGraph::generateGraph()
{
    while(kmer_list.size())
    {
        uint64_t kmer = kmer_list.back();
        kmer_list.pop_back();//destroy as we go along to save memory

        std::cout << kmer << std::endl;

        //generate prefix as all but the last 4 bits
        uint64_t prefix = (kmer & (PSFIX_MASK << 2)) >> 2;

        //generate suffix as all but the first base (2 bits per base)
        uint64_t suffix = kmer & PSFIX_MASK;

        bool prefix_revcmp,suffix_revcmp;
        EulerNode*prefix_node = nullptr;
        EulerNode*suffix_node = nullptr;

        //lookup or generate node for prefix and suffic sequences
        //convert to reverse complement if not already canonical
        generateNode(prefix,prefix_revcmp,prefix_node);//sets revcmp and node
        generateNode(suffix,suffix_revcmp,suffix_node);

        //create associated edge
        EulerEdge*edge = new EulerEdge(prefix_node,suffix_node,prefix_revcmp,suffix_revcmp);

        prefix_node->addOutgoing(edge);
        suffix_node->addIncoming(edge);
        graph_edges.push_back(edge);
    }
}

EulerEdge::EulerEdge(EulerNode*pfix_node,EulerNode*sfix_node,bool pfix_revcmp,bool sfix_revcmp)
:prefix_node(pfix_node),suffix_node(sfix_node),prefix_revcmp(pfix_revcmp),suffix_revcmp(sfix_revcmp)
{}

EulerNode::EulerNode(uint64_t psfix)
:sequence(psfix)
{
}

void EulerNode::addOutgoing(EulerEdge*edge)
{
    outgoing.push_back(edge);
}

void EulerNode::addIncoming(EulerEdge*edge)
{
    incoming.push_back(edge);
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
                rkmer =  (rkmer >> 2) + (uint64_t(3) << KMER_SHIFT);
                break;
            case 'C':
            case 'c':
                fkmer = ((fkmer << 2) +  uint64_t(1)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(2) << KMER_SHIFT);
                break;
            case 'G':
            case 'g':
                fkmer = ((fkmer << 2) +  uint64_t(2)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(1) << KMER_SHIFT);
                break;
            case 'T':
            case 't':
                fkmer = ((fkmer << 2) +  uint64_t(3)) & KMER_MASK;
                rkmer =  (rkmer >> 2) + (uint64_t(0) << KMER_SHIFT);
                break;
            default:
                throw std::runtime_error("Error: invalid character in kmer: " + seq[p]);
        }
    }

    //return canonical version
    if(fkmer <= rkmer) return fkmer;
    else                return rkmer;
}

//ensure psfix (prefix or suffix) is canonical, return true if it needed reverse complementing
bool psfix2canonical(uint64_t*psfix)
{
    uint64_t fwd_fix=*psfix,rev_fix=0;

    for(auto p=0; p<PSFIX_SIZE; p++)
    {
        //rev_fix: grows from the left (MSB)
        //complement the two least significant bits, shift then to the start
        rev_fix = (rev_fix >> 2) + ((~fwd_fix & 0x3) << PSFIX_SHIFT);

        //consume the used bits
        fwd_fix >>= 2;
    }

    //if *psfix is already canonical nothing to do
    if(*psfix <= rev_fix) return false;

    //replace *psfix with the canonical form
    *psfix = rev_fix;

    return true;
}

} //namespace kmerz
