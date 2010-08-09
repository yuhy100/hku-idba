#ifndef __SEQUENCE_H_

#define __SEQUENCE_H_

#include "globals.h"
#include "Kmer.h"

#include <algorithm>
#include <cstring>
#include <string>
#include <iostream>

class Read;

class Sequence
{
public:
    friend std::ostream &operator <<(std::ostream &stream, const Sequence &seq);

    Sequence(): bases(0), length(0), capacity(0) {}

    Sequence(const Sequence &seq): bases(0), length(0), capacity(0) 
    { SetContent(seq.bases, seq.length); }

    Sequence(const std::string &s): bases(0), length(0), capacity(0) 
    { SetContent(s.c_str(), s.size()); }

    Sequence(const char *s): bases(0), length(0), capacity(0)
    { SetContent(s, std::strlen(s)); }

    Sequence(const Read &read): bases(0), length(0), capacity(0) 
    { SetContent(read); }

    Sequence(const Kmer &kmer): bases(0), length(0), capacity(0) 
    { SetContent(kmer, kmerLength); }

    ~Sequence() { delete [] bases; }
    void Clear() { length = 0; }

    void Swap(Sequence &seq)
    {
        if (this != &seq)
        {
            std::swap(bases, seq.bases);
            std::swap(length, seq.length);
            std::swap(capacity, seq.capacity);
        }
    }

    char &operator [](int index) { return bases[index]; }
    const char &operator [](int index) const { return bases[index]; }

    const Sequence &operator =(const Sequence &seq)
    {
        if (this != &seq) 
            SetContent(seq.bases, seq.length);
        return *this;
    }

    const Sequence &operator =(const std::string &s)
    { SetContent(s.c_str(), s.size()); return *this; }

    const Sequence &operator =(const char *seq)
    { SetContent(seq, std::strlen(seq)); return *this; }

    const Sequence &operator =(const Read &read)
    { SetContent(read); return *this; }

    const Sequence &operator =(const Kmer &kmer)
    { SetContent(kmer, kmerLength); return *this; }

    const Sequence &operator +=(const Sequence &seq)
    { AddSequence(seq.bases, seq.length); return *this; }

    const Sequence &operator +=(const std::string &s)
    { AddSequence(s.c_str(), s.size()); return *this; }

    const Sequence &operator +=(const char *s)
    { AddSequence(s, std::strlen(s)); return *this; }

    const Sequence &operator +=(unsigned value)
    { AddNucleotide(value); return *this; }

    bool operator ==(const Sequence &seq) const
    {
        if (length != seq.length)
            return false;
        return std::memcmp(bases, seq.bases, length) == 0;
    }

    bool operator !=(const Sequence &seq) const
    {
        if (length != seq.length)
            return true;
        return std::memcmp(bases, seq.bases, length) != 0;
    }

    bool operator <(const Sequence &seq) const
    {
        int len = std::min(length, seq.length);
        int diff = memcmp(bases, seq.bases, len);
        if (diff != 0)
            return diff < 0;
        else
            return length < seq.length;
    }

    char *ToCString() { return bases; }
    const char *ToCString() const { return bases; }

    bool IsChar() const;
    bool IsCodon() const;

    int Size() const { return length; }
    int Capacity() const { return capacity; }
    void Resize(int l) { length = l; bases[l] = '\0'; }

    void Encode();
    void Decode();
    void ReverseComplement();

    void Trim(int t) { length -= t; bases[length] = '\0'; }
    void TrimError();

    void GetSubSequence(Sequence &seq, int offset, int subLength) const;
    Kmer GetKmer(int offset, int kmerLength) const;

protected:
    void SetContent(const Sequence &seq) { SetContent(seq.bases, seq.length); }
    void SetContent(const std::string &s) { SetContent(s.c_str(), s.size()); }
    void SetContent(const char *seq, int newLength);
    void SetContent(const Read &read);
    void SetContent(const Kmer &kmer, int kmerLength);

    void AddSequence(const Sequence &seq)
    { AddSequence(seq.bases, seq.length); }

    void AddSequence(const std::string &s)
    { AddSequence(s.c_str(), s.size()); }

    void AddSequence(const char *seq, int length)
    {
        if (capacity < this->length + length + 1)
        {
            Reallocate(this->length + length + 1);
        }
        std::copy(seq, seq + length, bases + this->length);
        this->length += length;
        bases[this->length] = '\0';
    }

    void AddNucleotide(unsigned value)
    {
        if (capacity == 0)
            Reallocate(4);
        else if (length+2 > capacity)
            Reallocate(capacity<<1);
        bases[length++] = value;
        bases[length] = '\0';
    }

private:
    void Reallocate(int new_capacity);

    static int code[256];

    char *bases;
    int length;
    int capacity;
};

#endif

