#ifndef __SEQUENCE_H_

#define __SEQUENCE_H_

#include "globals.h"
#include "kmer.h"
#include "read.h"

#include <cstdio>
#include <algorithm>

struct CompactSequence;
struct Read;

struct Sequence
{
    char *bases;
    int length;
    int capacity;

    Sequence() { bases = NULL; length = 0; capacity = 0; }
    Sequence(char *seq, int length)
    {
        bases = 0;
        length = 0;
        capacity = 0;
        SetContent(seq, length);
    }
    Sequence(const Sequence &seq)
    {
        bases = 0;
        length = 0;
        capacity = 0;
        SetContent(seq.bases, seq.length);
    }
    ~Sequence() { delete [] bases; }

    operator char *() { return bases; }
    char &operator [](int index) { return bases[index]; }
    const Sequence &operator =(const Sequence &seq)
    {
        SetContent(seq.bases, seq.length);
        return *this;
    }
    const Sequence &operator +=(const Sequence &seq)
    {
        if (capacity < length + seq.length + 1)
        {
            Reallocate(capacity < length + seq.length + 1);
        }
        std::copy(seq.bases, seq.bases + seq.length, bases + length);
        length += seq.length;
        bases[length] = '\0';
        return *this;
    }
    const Sequence &operator +=(unsigned value)
    {
        if (capacity == 0)
            Reallocate(4);
        else if (length+2 > capacity)
            Reallocate(capacity<<1);
        bases[length++] = value;
        bases[length] = '\0';
        return *this;
    }

    void Clear() { length = 0; }
    void Reallocate(int newCapacity);
    void SetContent(char *seq, int newLength);
    void SetContent(const Kmer &kmer, int kmerLength);
    void SetContent(CompactSequence *seq);
    void SetContent(Read *read);
    void SetContent(Sequence *seq) { SetContent(seq->bases, seq->length); }

    void AddNucleotide(unsigned value)
    {
        if (capacity == 0)
            Reallocate(4);
        else if (length+2 > capacity)
            Reallocate(capacity<<1);
        bases[length++] = value;
        bases[length] = '\0';
    }

    bool IsChar();
    bool IsCodon();
    void Encode();
    void Decode();

    void Resize(int l) { length = l; bases[l] = '\0'; }
    void Trim(int t) { length -= t; bases[length] = '\0'; }

    void ReverseComplement();
    void GetSubSequence(Sequence *seq, int offset, int subLength);
    Kmer GetKmer(int offset, int kmerLength);
};

#endif
