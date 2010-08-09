#ifndef __KMER_VECTOR_H_

#define __KMER_VECTOR_H_

#include "globals.h"
#include "Sequence.h"
#include "Kmer.h"

#include <algorithm>
#include <vector>

using namespace std;

class KmerVector
{
public:
    KmerVector(int length = 4)
    {
        this->length = length;
        this->size = (1 << (length<<1))/2 + (~length&1) * (1 << length)/2;
        values.resize(size);
        ranks.resize(size);
    }
    ~KmerVector() { }

    void Clear()
    {
        for (int i = 0; i < size; ++i)
            values[i] = 0;
    }

    KmerVector(const KmerVector &x)
    {
        length = x.length;
        size = x.size;
        values = x.values;
        ranks = x.ranks;
    }

    const KmerVector &operator =(const KmerVector &x)
    {
        length = x.length;
        size = x.size;
        values = x.values;
        ranks = x.ranks;
        return *this;
    }

    const KmerVector &operator +=(const KmerVector &x)
    {
        for (int i = 0; i < size; ++i)
            values[i] += x.values[i];
        ComputeRank();

        return *this;
    }

    const KmerVector &operator /=(double x)
    {
        for (int i = 0; i < size; ++i)
            values[i] /= x;

        return *this;
    }

    void ComputeRank();
    void Compute(Sequence *seq);

    int length;
    int size;
    std::vector<double> values;
    std::vector<int> ranks;
};

#endif
