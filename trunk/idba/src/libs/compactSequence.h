#ifndef __COMPACTSEQUENCE_H_

#define __COMPACTSEQUENCE_H_

#include "globals.h"

struct Sequence;

struct CompactSequence
{
    unsigned char *compressed;
    int length;
    int capacity;

    CompactSequence() { compressed = 0; length = 0; capacity = 0; }
    ~CompactSequence() { delete [] compressed; }

    void Clear() { length = 0; }
    void Reallocate(int newCapcity);
    void SetContent(Sequence *seq);

    unsigned GetNucleotide(unsigned index)
    {
        return (compressed[index>>2] >> ((index&3) << 1)) & 3;
    }
    void SetNucleotide(unsigned index, unsigned value)
    {
        compressed[index>>2] &= (unsigned char)~(3 << ((index&3) << 1));
        compressed[index>>2] |= value << ((index&3) << 1);
    }
    void AddNucleotide(unsigned value)
    {
        if ((((length+1) >> 2) + 1) > capacity)
        {
            Reallocate(capacity << 1);
        }
        SetNucleotide(length++, value);
    }

    void ReverseComplement();
};

#endif
