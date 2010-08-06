#ifndef __HASHGRAPH_H_

#define __HASHGRAPH_H_

#include "globals.h"
#include "hashNode.h"

#include <algorithm>
#include <iostream>
#include <vector>

const unsigned HashGraphFlagDeadend =   0x00001000;
const unsigned HashGraphFlagAlive =     0x00002000;
const unsigned HashGraphFlagUsed =      0x00004000;
const unsigned HashGraphFlagScaning =   0x00000001;

struct Sequence;
struct Kmer;

struct HashGraph
{
//     unsigned char *sync;
    HashNode **table;
    uint64 tableSize;
    uint64 numNodes;
    uint64 numEdges;
    unsigned minLength;

    Sequence corrected;
    int correctedWeight;

    HashGraph(uint64 tableSize);
    ~HashGraph();

    void Reallocate(uint64 newTableSize);

    void Clear();
    void ClearGraph();

    HashNode *InsertKmer(const Kmer &kmer, int count = 1);
    HashNode *GetByKmer(const Kmer &kmer);

    void InsertSequence(const Sequence *seq);
    void Chop(Sequence *seq);
    bool IsValid(Sequence *seq);
    bool IsContainBranch(Sequence *seq);
    bool CheckEdges(Sequence *seq);

    int GetDataFromSequence(const Sequence *seq);
    int AddDataFromSequence(const Sequence *seq, int data);
    int AddEdgesFromSequence(const Sequence *seq);
    int AddEdgesFromSequence(const Sequence *seq, unsigned char *exist);
    int AddEdgesFromContig(const Sequence *contig);

    uint64 RemoveDeadEnd(unsigned minCount, unsigned from = 0, unsigned to = 0, unsigned scan = HashGraphFlagScaning);
    int CheckDeadEnd(const Kmer &kmer, unsigned remain,
                     unsigned scan = HashGraphFlagScaning);
    void UpdateDepth(const Kmer &kmer);

    void Refresh(unsigned minCount = 1);
    void RefreshVertices(unsigned from = 0, unsigned to = 0, unsigned minCount = 1);
    void RefreshEdges(unsigned from = 0, unsigned to = 0);
    bool ErrorCorrect(Sequence *seq, int index, Kmer kmer, unsigned error);
    void ErrorCorrectWeight(Sequence *seq, int index, Kmer kmer, unsigned error, int sum);

//     bool FindPath(Sequence *path, const Kmer &from, const Kmer &to, int step);

    uint64 RemoveBubble();

    double AverageCoverage();
//     uint64 RemoveLowCoverageContigs();
    uint64 RemoveLowCoverageContigs(double c);
    uint64 Assemble(Sequence *contigs);
    uint64 Assemble(std::vector<Sequence> &contigs,
                    std::vector<Kmer> &branches);
};

#endif
