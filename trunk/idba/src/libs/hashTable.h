#ifndef __HASHGRAPH_H_

#define __HASHGRAPH_H_

#include "globals.h"
#include "hashNode.h"

struct Sequence;
struct Kmer;

struct HashTable
{
//     unsigned char *sync;
    HashNode **table;
    uint64 tableSize;
    uint64 numNodes;
    bool isReversable;

    HashTable(uint64 tableSize, bool isReversable = false);
    ~HashTable();

    void Clear();
    HashNode *InsertKmer(const Kmer &kmer);
    HashNode *GetByKmer(const Kmer &kmer);
};

#endif
