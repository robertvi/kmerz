#ifndef __ROBERTVI_KMERZ_CONFIG_H__
#define __ROBERTVI_KMERZ_CONFIG_H__

#include <string>

namespace kmerz
{
    class KmerzConfig
    {
        public:
            KmerzConfig(int argc,char**argv);
            virtual ~KmerzConfig();
        private:
            int minKmerCount;
            std::string inputFile;
            std::string outputFile;
    };
}

#endif //__ROBERTVI_KMERZ_PARSEARGS_H__
