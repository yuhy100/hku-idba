#include "globals.h"
#include "Log.h"
#include "Sequence.h"
#include "Utils.h"
#include "HashGraph.h"
#include "ContigGraph.h"
#include "HashAlign.h"
#include "MapGraph.h"
#include "Reader.h"

#include <boost/program_options.hpp>

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>

namespace po = boost::program_options;

using namespace std;

HashGraph hashGraph1;
HashGraph hashGraph2;
char line[MaxLine];
string comment1;
string comment2;
char buf[MaxLine];

void Usage()
{
    printf("Usage: compareKmer ref-file1 ref-file2\n");
}

int64 CommonKmer(HashGraph *graph1, HashGraph *graph2)
{
    int64 common = 0;
    for (HashGraph::Iterator iter = graph1->Begin(); iter != graph1->End(); ++iter)
    {
        KmerNode *node = *iter;
        if (graph2->GetNode(node->GetKmer()) != NULL)
            ++common;
    }
    return common;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        Usage();
        exit(1);
    }

    AddParameter("kmer", &kmerLength, INTEGER);
    ProcessParameters(argc, argv);

    FastAReader reader1(argv[1]);
    FastAReader reader2(argv[2]);

    Sequence ref1, ref2;

    reader1.Read(ref1, comment1);
    reader2.Read(ref2, comment2);

    ref1.Encode();
    ref2.Encode();

    hashGraph1.InsertSequence(ref1);

    while (reader1.Read(ref1, comment1))
        hashGraph1.InsertSequence(ref1);

    hashGraph2.InsertSequence(ref2);

    while (reader2.Read(ref2, comment2))
        hashGraph2.InsertSequence(ref2);

    int64 common = CommonKmer(&hashGraph1, &hashGraph2);

    printf("%.2f%%\n", 100.0 * 2 * common / (hashGraph1.NumNodes() + hashGraph2.NumNodes()));

//    printf("kmer: %d\n", kmerLength);
//
//    printf("ref1: %lld\n", hashGraph1.NumNodes());
//    printf("ref2: %lld\n", hashGraph2.NumNodes());
//    int64 common = CommonKmer(&hashGraph1, &hashGraph2);
//    printf("common: %lld only in ref1: %lld only in ref2: %lld\n", 
//            common, hashGraph1.NumNodes() - common, hashGraph2.NumNodes() - common);

    return 0;
}
