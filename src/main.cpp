#include "config.h"
#include "kmerz.h"

#include <iostream>

int main(int argc,char**argv)
{
    //parse command line arguments
    kmerz::KmerzConfig config(argc,argv);

    //load and min-count filter kmers from file
    kmerz::EulerGraph euler( config.getInputFile(),config.getMinKmerCount() );

    std::vector< std::string > contig_list;

    //generate contigs by walking the (implicit) overlap graph
    euler.generateContigs( contig_list );

    for(int i=0; i<contig_list.size(); i++)
    {
        std::cout << ">contig" + std::to_string(i) << std::endl;
        std::cout << contig_list.at(i) << std::endl;
    }

    return 0;
}
