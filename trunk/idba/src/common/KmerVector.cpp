#include "globals.h"
#include "Kmer.h"
#include "KmerVector.h"

#include <iostream>
#include <algorithm>
#include <cstdio>

using namespace std;

static long long table[1<<20];
static int aux[MaxLine];

class Comparer
{
public:
    Comparer(KmerVector *kv) { this->kv = kv; }

    bool operator ()(int x, int y) const 
    { return kv->values[x] < kv->values[y]; }

private:
    KmerVector *kv;
};

void KmerVector::ComputeRank()
{
    for (int i = 0; i < size; ++i)
        aux[i] = i;

    sort(aux, aux + size, Comparer(this));
    for (int i = 0; i < size; ++i)
        ranks[aux[i]] = i;
}

void KmerVector::Compute(Sequence *seq)
{
    int tableSize = 1 << (length << 1);
    fill_n(table, tableSize, 0);

    uint64 kmer = 0;
    for (int i = 0; i < length-1; ++i)
    {
        kmer |= (uint64((*seq)[i]) << (i << 1));
    }

    // Count the kmers
    for (int i = length-1; i < seq->Size(); ++i)
    {
        kmer |= (uint64((*seq)[i]) << ((length-1) << 1));
        ++table[kmer];
        kmer >>= 2;
    }

    // Sum up the count of the kmer and its reverseComplement,
    // and build the kmer vector.
    int k = 0;
    for (int i = 0; i < tableSize; ++i)
    {
        uint64 j = i;
        BitOperation::ReverseComplement(j);
        j >>= (64 - (length << 1));
        if (table[i] >= 0)
            values[k++] = (table[i] + table[j])*1.0/seq->Size();
        table[i] = table[j] = -1;
    }

    ComputeRank();
}
