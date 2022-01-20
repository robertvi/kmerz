#ifndef __ROBERTVI_KMERZ_CONFIG_H__
#define __ROBERTVI_KMERZ_CONFIG_H__

#include <string>

namespace kmerz
{
    class KmerzConfig
    {
        public:
            KmerzConfig(int _argc,char**_argv);
            virtual ~KmerzConfig();
            std::string getInputFile() const {return inputFile;}
            std::string getOutputFile() const {return outputFile;}
            int getMinKmerCount() const {return minKmerCount;}
        private:
            int findOption(const std::string&option,const std::string & alt="");
            std::string getOption(const std::string&option,const std::string & alt="");
            void parseArgs();
            void showOptions();

            //command line options
            int minKmerCount;
            std::string inputFile;
            std::string outputFile;

            int argc;
            char**argv;

            static const char*usage;
    };
} //namespace kmerz

#endif //__ROBERTVI_KMERZ_PARSEARGS_H__
