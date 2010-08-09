#ifndef __ASSEMBLY_ALGORITHM_H_

#define __ASSEMBLY_ALGORITHM_H_

#include "globals.h"
#include "Kmer.h"
#include "HashGraph.h"
#include "Read.h"
#include "Sequence.h"
#include "Contig.h"
#include "MapGraph.h"

#include <string>
#include <vector>

class AssemblyAlgorithm
{
public:
    AssemblyAlgorithm() {}

    int64 LoadReads(const std::string &readfile);

    std::vector<Sequence> longReads;
    std::vector<Read> reads;
    int trim;

private:
};

#endif
