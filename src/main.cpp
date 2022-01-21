#include "config.h"
#include "kmerz.h"

int main(int argc,char**argv)
{
    //parse command line arguments
    kmerz::KmerzConfig config(argc,argv);

    //load and min-count filter kmers from file
    kmerz::EulerGraph euler(config.getInputFile(),config.getMinKmerCount());

    //generate graph from kmer list
    euler.generateGraph();

    //generate assembled sequence(s) by walking through the graph
    euler.generatePaths(config.getOutputFile());

    return 0;
}
