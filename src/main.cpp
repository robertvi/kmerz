#include "config.h"
#include "kmerz.h"

int main(int argc,char**argv)
{
    kmerz::KmerzConfig config(argc,argv);

    kmerz::EulerGraph euler(config.getInputFile(),config.getMinKmerCount());

    return 0;
}
