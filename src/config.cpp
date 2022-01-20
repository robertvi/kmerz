#include "config.h"

#include <iostream>

namespace kmerz
{

KmerzConfig::KmerzConfig(int argc,char**argv)
{
    inputFile = "";
    outputFile = "-";

    std::cout << outputFile << std::endl;
}

KmerzConfig::~KmerzConfig()
{
}

}
