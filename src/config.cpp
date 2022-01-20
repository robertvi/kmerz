#include "config.h"

#include <iostream>
#include <stdexcept>

namespace kmerz
{

const char*KmerzConfig::usage =
"kmerz [--output|-o <filename>] [--min-count|-m <minKmerCount>] --input|-i <inputfile>";

KmerzConfig::KmerzConfig(int _argc,char**_argv)
:argc(_argc),argv(_argv)
{
    //set default option values
    inputFile = "";   //no default, must specify
    outputFile = "-"; //stdout
    minKmerCount = 0; //no filtering

    if(findOption("--help","-h"))
    {
        std::cout << "Usage: " << usage << std::endl;
        exit(0);
    }

    if(findOption("--output","-o"))
    {
        outputFile = getOption("--output","-o");
    }

    if(findOption("--min-count","-m"))
    {
        minKmerCount = std::stoi( getOption("--min-count","-m") );
    }

    if(findOption("--input","-i"))
    {
        inputFile = getOption("--input","-i");
    }
    else
    {
        throw std::runtime_error("Error: --input option is mandatory");
    }

    std::cout << "outputFile=" << outputFile << std::endl;
    std::cout << "inputFile="  << inputFile  << std::endl;
    std::cout << "minKmerCount=" << minKmerCount << std::endl;
}

KmerzConfig::~KmerzConfig()
{
}

//return position of first occurrence the option or 0 if not present
int KmerzConfig::findOption(const std::string&option,const std::string&alt)
{
    for(int i=1; i<argc; i++)
    {
        if(option.compare(argv[i]) == 0 || alt.compare(argv[i]) == 0) return i;
    }

    return 0;
}

//return string of the argv following the requested option or throw exception
std::string KmerzConfig::getOption(const std::string&option,const std::string&alt)
{
    for(int i=1; i<argc; i++)
    {
        if(option.compare(argv[i]) == 0 || alt.compare(argv[i]) == 0)
        {
            if(i < argc-1)
            {
                return argv[i+1];
            }

            throw std::runtime_error("Error: " + option + " has no associated value");
        }
    }

    throw std::runtime_error("Error: " + option + " not found");
}

} //kmerz namespace
