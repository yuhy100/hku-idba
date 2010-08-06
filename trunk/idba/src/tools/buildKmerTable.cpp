#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "utils.h"
#include "hashGraph.h"

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>

using namespace std;

vector<Sequence> longReads;
Read *reads;
char line[MaxLine];
char comment[MaxLine];
unsigned prefixLength = 3;
unsigned minCount = 2;
// unsigned trim = 0;
uint64 mask;
uint64 numReads;

void InsertSequence(HashGraph *hashGraph, const Sequence *seq, uint64 prefix)
{
    Kmer kmer = {{0}};
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(seq->bases[i]);
    for (int i = kmerLength-1; i < seq->length; ++i)
    {
        kmer.AddRight(seq->bases[i]);
        Kmer key = kmer;
        Kmer revComp = kmer;
        revComp.ReverseComplement();

        unsigned left, right;
        left = 4;
        right = 4;
        if (i >= (int)kmerLength)
            left = 3 - seq->bases[i-kmerLength];
        if (i+1 < seq->length)
            right = seq->bases[i+1];

        if (revComp < kmer)
        {
            key = revComp;
            swap(left, right);
        }

        if ((key.Hash() & mask) == prefix)
        {
            HashNode *node = hashGraph->InsertKmer(key);
            if (left != 4)
                node->in |= 1 << left;
            if (right != 4)
                node->out |= 1 << right;
        }
    }
}

int main(int argc, char *argv[])
{
    AddParameter("kmer", &kmerLength, INTEGER);
    AddParameter("prefixLength", &prefixLength, INTEGER);
    AddParameter("minCount", &minCount, INTEGER);
//     AddParameter("trim", &trim, INTEGER);

    ProcessParameters(argc, argv);

    LogMessage("Kmer %d\n", kmerLength);
    HashGraph hashGraph(MaxHashTable);

    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: buildKmerTable read-file table-file\n");
        fprintf(stderr, "       [--kmer k] [--minCount m] --prefixLength p]\n");
        exit(1);
    }

    mask = (1 << prefixLength) - 1;
    FILE *freadFile = OpenFile(argv[1], "rb");
    FILE *ftableFile = OpenFile(argv[2], "wb");

    Sequence seq;
    numReads = CountSequences(freadFile);
    reads = new Read[numReads];

    uint64 index = 0;
    while (ReadFasta(freadFile, &seq))
    {
        if (seq.IsChar())
        {
            seq.Encode();
//             seq.Trim(trim);

            if (seq.length < (int)kmerLength)
                continue;

            if (seq.length <= (int)MaxReadLength)
                reads[index++].SetContent(&seq);
            else
                longReads.push_back(seq);
        }
    }

    numReads = index;
    uint64 totalKmer = 0;
    uint64 maximum = 0;
    uint64 remain = 0;
    for (int prefix = 0; prefix < (1 << prefixLength); ++prefix)
    {
        hashGraph.Clear();
        for (uint64 i = 0; i < numReads; ++i)
        {
            seq.SetContent(&reads[i]);
            InsertSequence(&hashGraph, &seq, prefix);
        }

        for (unsigned i = 0; i < longReads.size(); ++i)
            InsertSequence(&hashGraph, &longReads[i], prefix);

        LogMessage("reads %d total %d\n", numReads, hashGraph.numNodes);
        remain += WriteHashGraph(ftableFile, &hashGraph, minCount);

        totalKmer += hashGraph.numNodes;
        if (hashGraph.numNodes > maximum)
            maximum = hashGraph.numNodes;
    }

    LogMessage("read: %s\n", argv[1]);
    LogMessage("max/tol %d/%d %.2f%%\n", maximum, totalKmer, 100.0*maximum/totalKmer);
    LogMessage("rem/tol %d/%d %.2f%%\n\n", remain, totalKmer, 100.0*remain/totalKmer);

    fclose(freadFile);
    fclose(ftableFile);

    return 0;
}
