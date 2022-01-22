#include "kmerz.h"

#include <iostream>
#include <fstream>

const uint64_t KMER_SIZE= 31;
const std::string bases = "ATCG";

namespace kmerz
{

EulerGraph::EulerGraph(std::vector< std::string >&input_list)
{
    for(auto it=input_list.begin(); it!=input_list.end(); it++)
    {
        std::string kmer = *it;
        makeCanonical( kmer );
        kmer_set.emplace( kmer );
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

char EulerGraph::extendSuffix(std::string suffix)
{
    char next_base = 'X';
    bool unique = true;
    std::string kmer;

    //find all possible extensions and remove them from the kmer dictionary
    for(int i=0; i<4; i++)
    {
        kmer.assign(suffix);
        kmer.append(1, char(bases[i]));
        makeCanonical(kmer);

        //search for the possible kmer extending the suffix
        auto set_it = kmer_set.find(kmer);

        //if not present continue to next base
        if(set_it == kmer_set.end()) continue;

        //remove the kmer from the set
        kmer_set.erase(set_it);

        if(next_base == 'X')
        {
            //possible unique extension
            next_base = char(bases[i]);
        }
        else
        {
            //non unique extension
            unique = false;
        }
    }

    if(next_base != 'X' && unique == true) return next_base;

    return 'X';
}

std::string EulerGraph::walkForwards(std::string kmer)
{
    std::string contig = "";

    while(true)
    {
        kmer.assign(kmer.substr(1));

        char next_base = extendSuffix(kmer);

        //stop if no valid extension found
        if(next_base == 'X') return contig;

        //extend contig and kmer with next base
        contig.append(1,next_base);
        kmer.append(1,next_base);
    }
}

void EulerGraph::generateContigs(std::vector< std::string >&contig_list)
{
    //while kmers left in set
    while(kmer_set.size())
    {
        //pick arbitrary seed kmer
        auto seed_it = kmer_set.begin();
        std::string seed = *seed_it;
        std::string rev_seed = reverseComplement(seed);

        //remove the seed kmer to mark it as visited already
        kmer_set.erase(seed_it);

        //walk sequence forward from seed
        std::string fwd_contig = walkForwards(seed);

        //walk sequence backward from seed
        std::string rev_contig = walkForwards(rev_seed);

        //join into final contig
        std::string contig = reverseComplement(rev_contig) + seed + fwd_contig;

        contig_list.push_back(contig);

        /*std::cout << "fwd:" << std::endl;
        std::cout << fwd_contig << std::endl;
        std::cout << "rev:" << std::endl;
        std::cout << rev_contig << std::endl;
        std::cout << "contig:" << std::endl;
        std::cout << contig << std::endl;*/
    }
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

std::string getSuffix(const std::string&seq)
{
    return seq.substr(1);
}

} //namespace kmerz
