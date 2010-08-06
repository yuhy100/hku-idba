#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "compactSequence.h"

#include <cstdio>
#include <algorithm>
#include <iostream>

using namespace std;

void CompactSequence::Reallocate(int newCapacity)
{
#ifdef DEBUG
    if (newCapacity < (length>>2) + 1)
    {
        LogError("CompactSequence::Reallocate() The new capacity is too short\n");
        exit(1);
    }
#endif
    unsigned char *newCompressed = new unsigned char[newCapacity];
    if (newCompressed == NULL)
    {
        LogError("CompactSequence::Reallocate() No enough memory\n");
        exit(1);
    }

    if (length != 0)
    {
        copy(compressed, compressed + (length>>2) + 1, newCompressed);
    }
    capacity = newCapacity;
    delete [] compressed;
    compressed = newCompressed;
}

void CompactSequence::SetContent(Sequence *seq)
{
    if (!seq->IsCodon())
    {
        LogError("CompactSequence::SetSequence: not in codon format\n");
        exit(1);
    }

    Clear();
    if ((seq->length>>2) + 1 > capacity
         || (seq->length>>2) + 1 < (capacity << 1))
        Reallocate((seq->length>>2) + 1);
    fill_n(compressed, capacity, 0);
    for (int i = 0; i < seq->length; ++i)
    {
        compressed[i>>2] |= seq->bases[i] << ((i&3) << 1);
    }
    length = seq->length;
}

void CompactSequence::ReverseComplement()
{
    for (int i = 0, j = length-1; i <= j; ++i, --j)
    {
        if (i != j)
        {
            int x = GetNucleotide(i);
            int y = GetNucleotide(j);
            SetNucleotide(i, 3 - y);
            SetNucleotide(j, 3 - x);
        }
        else
        {
            SetNucleotide(i, 3 - GetNucleotide(i));
        }
    }
}