#ifndef __SIMPLE_ALGORITHM_

#define __SIMPLE_ALGORITHM_

#include "globals.h"
#include "Kmer.h"
#include "HashGraph.h"
#include "Read.h"
#include "Sequence.h"
#include "Contig.h"
#include "MapGraph.h"
#include "AssemblyAlgorithm.h"
#include "Reader.h"

#include <string>
#include <vector>

class SimpleAlgorithm: public AssemblyAlgorithm
{
public:
    SimpleAlgorithm();

    int64 GenerateKmers();
    void AddbackInternalKmers();
    void RemoveLowCoverageContigs();
    void AddEdgesToContigGraph();
    void Iterate();

    void IDBA(FILE *fread, std::vector<Contig> &contigs);

    void Run();

private:
public:
    static const uint64 MaxDistance = 1000000;

    std::string prefix;
    std::string readfile;
    std::string contigfile;
    std::string scaffile;
    int mink;
    int maxk;
    int minLength;
    HashGraph hashGraph;
    unsigned prefixLength;
    int minCount;
    uint64 mask;
    int64 numReads;
    std::vector<Contig> contigs;
    std::vector<Kmer> branches;
    double cover;
    bool isScaffold;

    int totalContigs;
    int alignedReads;
    HashAlign hashAlign;
    std::vector<std::vector<Alignment> > alignResults;

    int minPairs;
    int hits[MaxDistance];
    int offset;
    MapGraph mapGraph;
    int estimateDistance;
    int estimateRange;
    int minContig;

};
#endif
