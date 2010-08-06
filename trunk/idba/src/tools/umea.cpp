#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "utils.h"
#include "hashGraph.h"

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;

vector<Sequence> contigs;
vector<int> length;
vector<Kmer> branches;
char line[MaxLine];
HashGraph *hashGraph = NULL;
unsigned minCount = 2;
unsigned minLength = 0;
unsigned minContig = 100;
int cover = 0;

int main(int argc, char *argv[])
{
    AddParameter("kmer", &kmerLength, INTEGER);
    AddParameter("minCount", &minCount, INTEGER);
//     AddParameter("minLength", &minLength, INTEGER);
    AddParameter("minContig", &minContig, INTEGER);
    AddParameter("cover", &cover, INTEGER);

    ProcessParameters(argc, argv);

    if (minLength == 0)
        minLength = kmerLength * 2;
    if (minLength > 255)
        minLength = 240;

    if (argc < 4
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: umea read-file table-file contig-file\n");
        fprintf(stderr, "       [--Kmer k] [--minCount m]\n");
        fprintf(stderr, "       [--minContig l] [--cover c] \n");
        exit(1);
    }

//     FILE *freadFile = OpenFile(argv[1], "rb");
    FILE *ftableFile = OpenFile(argv[2], "rb");
    FILE *fcontigFile = OpenFile(argv[3], "wb");

    hashGraph = new HashGraph(MaxHashTable);

    ReadHashGraph(ftableFile, hashGraph, minCount);
    fclose(ftableFile);

    Sequence seq;
    unsigned index = 0;
    hashGraph->RefreshVertices();
    hashGraph->RefreshEdges();

    LogMessage("total kmer: %d\t edges: %d\n",
        hashGraph->numNodes, hashGraph->numEdges);

    unsigned deadEnd = hashGraph->RemoveDeadEnd(minLength);
    LogMessage("Remove %d dead ends from %d nodes\n", deadEnd, hashGraph->numNodes);

    unsigned remain = hashGraph->numNodes - deadEnd;
    LogMessage("remain %d kmers\n", remain);

    hashGraph->RefreshVertices();
    hashGraph->RefreshEdges();

    int lowContigs = 0;
    if (cover == 0)
        lowContigs = hashGraph->RemoveLowCoverageContigs(hashGraph->AverageCoverage()/5);
    else
        lowContigs = hashGraph->RemoveLowCoverageContigs(cover);
    LogMessage("remove %d low coverage contigs\n", lowContigs);

    hashGraph->RefreshVertices();
    hashGraph->RefreshEdges();

    unsigned bubble = hashGraph->RemoveBubble();
    LogMessage("Remove %d bubbles\n", bubble);

    unsigned total = hashGraph->Assemble(contigs, branches);
    LogMessage("get contigs %d\n", total);

    index = 0;
    int contigs100 = 0;
    for (unsigned i = 0; i < total; ++i)
    {
        {
            ++contigs100;
            fprintf(fcontigFile, ">contig%d_length_%d\n",
                    index++, contigs[i].length);
            WriteFasta(fcontigFile, &contigs[i]);
        }
    }

    length.resize(contigs.size());
    for (unsigned i = 0; i < contigs.size(); ++i)
        length[i] = contigs[i].length;
    sort(length.begin(), length.end());
    reverse(length.begin(), length.end());

    int n50 = 0;
    long long sum = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        sum += length[i];
        if (sum >= 0.5*remain && n50 == 0)
        {
            n50 = length[i];
        }
    }

    printf("read: %s\nkmer: %d minCount: %d minLength: %d\n",
           argv[1], kmerLength, minCount, minLength);
    printf("total: %d N50: %d ave: %lld nodes: %llu\n",
           contigs100, n50, sum/contigs100, hashGraph->numNodes);

    fclose(fcontigFile);

    return 0;
}
