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
    euler.generateContigs(contig_list);

    //print to stdout as fasta format
    euler.printToStdoutFlat(contig_list);

    //success
    return 0;
}
