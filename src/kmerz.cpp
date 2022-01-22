#include "kmerz.h"

#include <iostream>
#include <fstream>

const uint64_t KMER_SIZE= 31;

namespace kmerz
{

EulerGraph::EulerGraph(std::vector< std::string >&input_list)
{
    for(auto it=input_list.begin(); it!=input_list.end(); it++)
    {
        std::string canonical_kmer = *it;
        makeCanonical( canonical_kmer );
        kmer_set.emplace( canonical_kmer );
    }
}

//load kmer sequences and counts from file
//format is kmer sequence, reverse complement of kmer sequence (ignored) and kmer count
//separated by spaces
//only the first of the two kmer sequences is used, converted to canonical if not already
//count is only used for simple min-mount filtering then discarded
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
        if(fwd_kmer.length() != KMER_SIZE)
        {
            throw std::runtime_error("Error: forward kmer not "+std::to_string(KMER_SIZE)+" bases on line: " + std::to_string(line_ctr) + ":"+fwd_kmer);
        }

        ifs >> rev_kmer;
        if(ifs.fail())
        {
            throw std::runtime_error("Error: failed to read reverse kmer on line: " + std::to_string(line_ctr));
        }
        if(rev_kmer.length() != KMER_SIZE)
        {
            throw std::runtime_error("Error: reverse kmer not "+std::to_string(KMER_SIZE)+" bases on line: " + std::to_string(line_ctr) + ":"+fwd_kmer);
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

        //filter out low coverage kmers
        if(count < minKmerCount) continue;

        //fwd_kmer should be canonical already
        if(makeCanonical(fwd_kmer))
        {
            throw std::runtime_error("Error: fwd_kmer on line: " + std::to_string(line_ctr) + " was not canonical");
        }

        kmer_set.emplace( fwd_kmer );
    }

    ifs.close();
}

void EulerGraph::generateContigs()
{

}

//make a string canonical
//return true if the value changed
bool makeCanonical(std::string&seq)
{
    std::string revcmp = reverseComplement(seq);

    if(revcmp.compare(seq) < 0)
    {
        seq.assign(revcmp);
        return true;
    }

    return false;
}

//reverse complement a string
std::string reverseComplement(const std::string&seq)
{
    std::string revcmp = "";

    for(int i=seq.size()-1; i>=0; i--)
    {
        switch(seq.at(i))
        {
            case 'A':
            case 'a':
                revcmp.append(1 ,'T');
                break;
            case 'T':
            case 't':
                revcmp.append(1 ,'A');
                break;
            case 'C':
            case 'c':
                revcmp.append(1 ,'G');
                break;
            case 'G':
            case 'g':
                revcmp.append(1 ,'C');
                break;
        }
    }

    return revcmp;
}


} //namespace kmerz
