#ifndef __KMER_H_

#define __KMER_H_

#include "globals.h"

#include <algorithm>
#include <iostream>

const uint64 SwapMask32 = 0x00000000FFFFFFFFULL;
const uint64 SwapMask16 = 0x0000FFFF0000FFFFULL;
const uint64 SwapMask8 = 0x00FF00FF00FF00FFULL;
const uint64 SwapMask4 = 0x0F0F0F0F0F0F0F0FULL;
const uint64 SwapMask2 = 0x3333333333333333ULL;

inline void ReverseComplement(uint64 &value)
{
    value = ((value & SwapMask32) << 32) | ((value & ~SwapMask32) >> 32);
    value = ((value & SwapMask16) << 16) | ((value & ~SwapMask16) >> 16);
    value = ((value & SwapMask8) << 8) | ((value & ~SwapMask8) >> 8);
    value = ((value & SwapMask4) << 4) | ((value & ~SwapMask4) >> 4);
    value = ((value & SwapMask2) << 2) | ((value & ~SwapMask2) >> 2);
    value = ~value;
}

struct Kmer
{
    uint64 data[KMERSIZE];

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
        ::ReverseComplement(data[0]);
        data[0] = data[0] >> (64 - (kmerLength << 1));
#elif KMERSIZE == 2
        if (kmerLength <= 32)
        {
            ::ReverseComplement(data[0]);
            data[0] = (data[0] >> (64 - (kmerLength << 1)));
        }
        else
        {
            ::ReverseComplement(data[0]);
            ::ReverseComplement(data[1]);
            std::swap(data[0], data[1]);
            data[0] = (data[0] >> (128 - (kmerLength << 1)))
                    | (data[1] << ((kmerLength << 1) - 64));
            data[1] = data[1] >> (128 - (kmerLength << 1));
        }
#elif KMERSIZE == 3
        if (kmerLength <= 32)
        {
                ::ReverseComplement(data[0]);
                data[0] = (data[0] >> (64 - (kmerLength << 1)));
        }
        else if (kmerLength <= 64)
        {
                ::ReverseComplement(data[0]);
                ::ReverseComplement(data[1]);
                std::swap(data[0], data[1]);
                data[0] = (data[0] >> (128 - (kmerLength << 1)))
                        | (data[1] << ((kmerLength << 1) - 64));
                data[1] = data[1] >> (128 - (kmerLength << 1));
        }
        else
        {
                ::ReverseComplement(data[0]);
                ::ReverseComplement(data[1]);
                ::ReverseComplement(data[2]);
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
                ::ReverseComplement(data[0]);
                data[0] = (data[0] >> (64 - (kmerLength << 1)));
        }
        else if (kmerLength <= 64)
        {
                ::ReverseComplement(data[0]);
                ::ReverseComplement(data[1]);
                std::swap(data[0], data[1]);
                data[0] = (data[0] >> (128 - (kmerLength << 1)))
                        | (data[1] << ((kmerLength << 1) - 64));
                data[1] = data[1] >> (128 - (kmerLength << 1));
        }
        else
        {
                ::ReverseComplement(data[0]);
                ::ReverseComplement(data[1]);
                ::ReverseComplement(data[2]);
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
            ::ReverseComplement(data[i]);
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
//         data[0] = (data[0] >> 2) | (data[1] << 62);
//         data[1] >>= 2;
//         data[(kmerLength-1) >> 5] |= uint64(ch) << (((kmerLength-1) & 31) << 1);
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
        return data[0] ^ data[1];
#elif KMERSIZE == 3
        return data[0] ^ data[1] ^ data[2];
#else
        uint64 key = 0;
        for (int i = 0; i < KMERSIZE; ++i)
            key ^= data[i];
        return key;
#endif
    }
};

#endif
