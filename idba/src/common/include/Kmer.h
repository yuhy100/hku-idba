#ifndef __KMER_H_

#define __KMER_H_

#include "globals.h"
#include "BitOperation.h"

#include <algorithm>
#include <iostream>

class Kmer32
{
public:
    friend std::ostream &operator <<(std::ostream &stream, const Kmer32 &kmer);
    friend class HashTable;

    Kmer32()
    { data[0] = 0; }

    bool operator <(const Kmer32 &kmer) const
    { return data[0] < kmer.data[0]; }

    bool operator ==(const Kmer32 &kmer) const
    { return data[0] == kmer.data[0]; }

    bool operator !=(const Kmer32 &kmer) const
    { return data[0] != kmer.data[0]; }

    void ReverseComplement()
    {
        BitOperation::ReverseComplement(data[0]);
        data[0] = data[0] >> (64 - (kmerLength << 1));
    }

    void AddRight(unsigned ch)
    {
        data[0] = (data[0] >> 2) | (uint64(ch) << ((kmerLength-1) << 1));
    }

    unsigned GetBase(unsigned index) const
    {
        index <<= 1;
        return (data[0] >> index) & 3;
    }

    void SetBase(unsigned index, unsigned ch)
    {
        index <<= 1;
        data[0] = (data[0] & ~(3ULL << index)) | (uint64(ch) << index);
    }

    uint64 Hash() const
    {
        return data[0];
    }

    bool IsPalindrome() const
    {
        Kmer32 kmer = *this;
        kmer.ReverseComplement();
        return *this == kmer;
    }

private:
    uint64 data[1];
};

class Kmer
{
public:
    friend std::ostream &operator <<(std::ostream &stream, const Kmer &kmer);
    friend class HashTable;

    Kmer()
    {
#if KMERSIZE == 1
        data[0] = 0;
#elif KMERSIZE == 2
        data[0] = data[1] = 0;
#elif KMERSIZE == 3
        data[0] = data[1] = data[2] = 0;
#else
        for (int i = 0; i < KMERSIZE; ++i)
            data[i] = 0;
#endif
    }

    bool operator <(const Kmer &kmer) const
    {
#if KMERSIZE == 1
        return data[0] < kmer.data[0];
#elif KMERSIZE == 2
        if (data[1] != kmer.data[1])
            return data[1] < kmer.data[1];
        else
            return data[0] < kmer.data[0];
#elif KMERSIZE == 3
        if (data[2] != kmer.data[2])
            return data[2] < kmer.data[2];
        else if (data[1] != kmer.data[1])
            return data[1] < kmer.data[1];
        else
            return data[0] < kmer.data[0];
#else
        for (int i = KMERSIZE-1; i >= 0; --i)
        {
            if (data[i] != kmer.data[i])
                return data[i] < kmer.data[i];
        }
        return false;
#endif
    }

    bool operator ==(const Kmer &kmer) const
    {
#if KMERSIZE == 1
        return data[0] == kmer.data[0];
#elif KMERSIZE == 2
        return data[0] == kmer.data[0] && data[1] == kmer.data[1];
#elif KMERSIZE == 3
        return data[0] == kmer.data[0] && data[1] == kmer.data[1] && data[2] == kmer.data[2];
#else
        for (int i = 0; i < KMERSIZE; ++i)
        {
            if (data[i] != kmer.data[i])
                return false;
        }
        return true;
#endif
    }

    bool operator !=(const Kmer &kmer) const
    {
#if KMERSIZE == 1
        return data[0] != kmer.data[0];
#elif KMERSIZE == 2
        return data[0] != kmer.data[0] || data[1] != kmer.data[1];
#elif KMERSIZE == 3
        return data[0] != kmer.data[0] || data[1] != kmer.data[1] || data[2] != kmer.data[2];
#else
        for (int i = 0; i < KMERSIZE; ++i)
        {
            if (data[i] != kmer.data[i])
                return true;
        }
        return false;
#endif
    }

    void ReverseComplement()
    {
#if KMERSIZE == 1
        BitOperation::ReverseComplement(data[0]);
        data[0] = data[0] >> (64 - (kmerLength << 1));
#elif KMERSIZE == 2
        if (kmerLength <= 32)
        {
            BitOperation::ReverseComplement(data[0]);
            data[0] = (data[0] >> (64 - (kmerLength << 1)));
        }
        else
        {
            BitOperation::ReverseComplement(data[0]);
            BitOperation::ReverseComplement(data[1]);
            std::swap(data[0], data[1]);
            data[0] = (data[0] >> (128 - (kmerLength << 1)))
                    | (data[1] << ((kmerLength << 1) - 64));
            data[1] = data[1] >> (128 - (kmerLength << 1));
        }
#elif KMERSIZE == 3
        if (kmerLength <= 32)
        {
            BitOperation::ReverseComplement(data[0]);
            data[0] = (data[0] >> (64 - (kmerLength << 1)));
        }
        else if (kmerLength <= 64)
        {
            BitOperation::ReverseComplement(data[0]);
            BitOperation::ReverseComplement(data[1]);
            std::swap(data[0], data[1]);
            data[0] = (data[0] >> (128 - (kmerLength << 1)))
                    | (data[1] << ((kmerLength << 1) - 64));
            data[1] = data[1] >> (128 - (kmerLength << 1));
        }
        else
        {
            BitOperation::ReverseComplement(data[0]);
            BitOperation::ReverseComplement(data[1]);
            BitOperation::ReverseComplement(data[2]);
            std::swap(data[0], data[2]);

            data[0] = (data[0] >> (192 - (kmerLength << 1)))
                    | (data[1] << ((kmerLength << 1) - 128));
            data[1] = (data[1] >> (192 - (kmerLength << 1)))
                    | (data[2] << ((kmerLength << 1) - 128));
            data[2] = data[2] >> (192 - (kmerLength << 1));
        }
#elif KMERSIZE == 3
        if (kmerLength <= 32)
        {
            BitOperation::ReverseComplement(data[0]);
            data[0] = (data[0] >> (64 - (kmerLength << 1)));
        }
        else if (kmerLength <= 64)
        {
            BitOperation::ReverseComplement(data[0]);
            BitOperation::ReverseComplement(data[1]);
            std::swap(data[0], data[1]);
            data[0] = (data[0] >> (128 - (kmerLength << 1)))
                    | (data[1] << ((kmerLength << 1) - 64));
            data[1] = data[1] >> (128 - (kmerLength << 1));
        }
        else
        {
            BitOperation::ReverseComplement(data[0]);
            BitOperation::ReverseComplement(data[1]);
            BitOperation::ReverseComplement(data[2]);
            std::swap(data[0], data[2]);

            data[0] = (data[0] >> (192 - (kmerLength << 1)))
                    | (data[1] << ((kmerLength << 1) - 128));
            data[1] = (data[1] >> (192 - (kmerLength << 1)))
                    | (data[2] << ((kmerLength << 1) - 128));
            data[2] = data[2] >> (192 - (kmerLength << 1));
        }
#else
        int l = (kmerLength >> 5) + ((kmerLength & 31) != 0);
        for (int i = 0; i < l; ++i)
            BitOperation::ReverseComplement(data[i]);
        for (int i = 0; i < (l>>1); ++i)
            std::swap(data[i], data[l-1-i]);
        int offset = ((l << 5) - kmerLength) << 1;
        for (int i = 0; i+1 < l; ++i)
            data[i] = (data[i] >> offset) | (data[i+1] << (64 - offset));
        data[l-1] >>= offset;
#endif
    }

    void AddRight(unsigned ch)
    {
#if KMERSIZE == 1
        data[0] = (data[0] >> 2) | (uint64(ch) << ((kmerLength-1) << 1));
#elif KMERSIZE == 2
        if (kmerLength <= 32)
        {
            data[0] = (data[0] >> 2) | (uint64(ch) << ((kmerLength-1) << 1));
        }
        else
        {
            data[0] = (data[0] >> 2) | (data[1] << 62);
            data[1] = (data[1] >> 2) | (uint64(ch) << ((kmerLength-33) << 1));
        }
#elif KMERSIZE == 3
        if (kmerLength <= 32)
        {
            data[0] = (data[0] >> 2) | (uint64(ch) << ((kmerLength-1) << 1));
        }
        else if (kmerLength <= 64)
        {
            data[0] = (data[0] >> 2) | (data[1] << 62);
            data[1] = (data[1] >> 2) | (uint64(ch) << ((kmerLength-33) << 1));
        }
        else
        {
            data[0] = (data[0] >> 2) | (data[1] << 62);
            data[1] = (data[1] >> 2) | (data[2] << 62);
            data[2] = (data[2] >> 2) | (uint64(ch) << ((kmerLength-65) << 1));
        }
#else
        int l = (kmerLength >> 5) + ((kmerLength & 31) != 0);
        for (int i = 0; i+1 < l; ++i)
            data[i] = (data[i] >> 2) | (data[i+1] << 62);
        data[l-1] = (data[l-1] >> 2) | (uint64(ch) << (((kmerLength - 1) & 31)<< 1));
#endif
    }

    unsigned GetBase(unsigned index) const
    {
#if KMERSIZE == 1
        index <<= 1;
        return (data[0] >> index) & 3;
#elif KMERSIZE == 2
        if (index < 32)
        {
            index <<= 1;
            return (data[0] >> index) & 3;
        }
        else
        {
            index = (index - 32) << 1;
            return (data[1] >> index) & 3;
        }
#elif KMERSIZE == 3
        if (index < 32)
        {
            index <<= 1;
            return (data[0] >> index) & 3;
        }
        else if (index < 64)
        {
            index = (index - 32) << 1;
            return (data[1] >> index) & 3;
        }
        else
        {
            index = (index - 64) << 1;
            return (data[2] >> index) & 3;
        }
#else
        return (data[index>>5] >> ((index & 31) << 1)) & 3;
#endif
    }

    void SetBase(unsigned index, unsigned ch)
    {
#if KMERSIZE == 1
        index <<= 1;
        data[0] = (data[0] & ~(3ULL << index)) | (uint64(ch) << index);
#elif KMERSIZE == 2
        if (index < 32)
        {
            index <<= 1;
            data[0] = (data[0] & ~(3ULL << index)) | (uint64(ch) << index);
        }
        else
        {
            index = (index & 31) << 1;
            data[1] = (data[1] & ~(3ULL << index)) | (uint64(ch) << index);
        }
#elif KMERSIZE == 3
        if (index < 32)
        {
            index <<= 1;
            data[0] = (data[0] & ~(3ULL << index)) | (uint64(ch) << index);
        }
        else if (index < 64)
        {
            index = (index - 32) << 1;
            data[1] = (data[1] & ~(3ULL << index)) | (uint64(ch) << index);
        }
        else
        {
            index = (index - 64) << 1;
            data[2] = (data[2] & ~(3ULL << index)) | (uint64(ch) << index);
        }
#else
        int offset = (index & 31) << 1;
        data[index>>5] = (data[index>>5] & ~(3ULL << offset)) | (uint64(ch) << offset);
#endif
    }

    uint64 Hash() const
    {
#if KMERSIZE == 1
        return data[0];
#elif KMERSIZE == 2
        return ((data[0] ^ data[1]) * 1299709 + 104729) % 323780508946331ULL;
#elif KMERSIZE == 3
        return data[0] ^ data[1] ^ data[2];
#else
        uint64 key = 0;
        for (int i = 0; i < KMERSIZE; ++i)
            key ^= data[i];
        return key;
#endif
    }

    bool IsPalindrome() const
    {
        Kmer kmer = *this;
        kmer.ReverseComplement();
        return *this == kmer;
    }

private:
    uint64 data[KMERSIZE];
};

#endif
