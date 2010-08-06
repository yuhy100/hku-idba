#include "globals.h"
#include "sequence.h"
#include "kmer.h"
#include "hashNode.h"
#include "hashTable.h"
#include "log.h"

#include <cstdio>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <queue>

using namespace std;

typedef HashNode *HashNodePointer;

HashTable::HashTable(uint64 tableSize, bool isReversable)
{
//     sync = new unsigned char[tableSize];
    table = new HashNodePointer[tableSize];
    this->tableSize = tableSize;
    this->isReversable = isReversable;
    numNodes = 0;

    if (table == NULL)
    {
        LogError("HashTable::HashTable() not enough memory\n");
        exit(1);
    }

    for (unsigned i = 0; i < tableSize; ++i)
    {
        table[i] = NULL;
    }
}

HashTable::~HashTable()
{
    Clear();
    delete [] table;
}

void HashTable::Clear()
{
    numNodes = 0;
    for (unsigned i = 0; i < tableSize; ++i)
    {
        HashNode *node = table[i];
        while (node != NULL)
        {
            HashNode *p = node;
            node = node->next;
            HashNode::FreeNode(p);
        }
        table[i] = NULL;
    }
}

HashNode *HashTable::InsertKmer(const Kmer &kmer)
{
    Kmer key = kmer;

    if (isReversable)
    {
        Kmer revComp = kmer;
        revComp.ReverseComplement();
        if (revComp < kmer)
            key = revComp;
    }

    uint64 index = key.Hash() % tableSize;
    HashNode *p = table[index];

    while (p != NULL && p->kmer != key)
        p = p->next;

    if (p == NULL)
    {
        p = HashNode::NewNode();
        p->SetContent(key);
        p->next = table[index];
        table[index] = p;
        ++numNodes;
    }

    ++p->count;

    return p;
}

HashNode *HashTable::GetByKmer(const Kmer &kmer)
{
    Kmer key = kmer;

    if (isReversable)
    {
        Kmer revComp = kmer;
        revComp.ReverseComplement();
        if (revComp < kmer)
            key = revComp;
    }

    uint64 index = key.Hash() % tableSize;
    HashNode *p = table[index];

    while (p != NULL && p->kmer != key)
        p = p->next;

    return p;
}

