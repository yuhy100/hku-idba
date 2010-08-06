#include "globals.h"
#include "log.h"
#include "kmer.h"
#include "compactSequence.h"
#include "sequence.h"
#include "read.h"

#include <cstdio>
#include <cstring>
#include <cctype>
#include <algorithm>

using namespace std;

static int code[] =
{
    'A', 'C', 'G', 'T', 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
    64, 0, 66, 1, 68, 69, 70, 2, 72, 73, 74, 75, 76, 77, 78, 79,
    80, 81, 82, 83, 3, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
    96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
};

// static int codeRevCmp[] =
// {
//     3, 2, 1, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
//     16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
//     32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
//     48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
//     64, 'T', 66, 'G', 68, 69, 70, 'C', 72, 73, 74, 75, 76, 77, 78, 79,
//     80, 81, 82, 83, 'A', 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
//     96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
//     112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
//     128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
//     144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
//     160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
//     176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
//     192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
//     208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
//     224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
//     240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
// };

void Sequence::Reallocate(int newCapacity)
{
#ifdef DEBUG
    if (newCapacity < length+1)
    {
        LogError("Sequence::Reallocate() The new capacity is too short\n");
        exit(1);
    }
#endif

    char *newBases = new char[newCapacity];
    if (newBases == NULL)
    {
        LogError("Sequence::Reallocate() No enough memory\n");
        exit(1);
    }

    if (length != 0)
    {
        copy(bases, bases + length+1, newBases);
    }
    capacity = newCapacity;
    delete [] bases;
    bases = newBases;
}

void Sequence::SetContent(char *seq, int newLength)
{
    Clear();
    if (newLength+1 > capacity || newLength+1 < (capacity>>1))
        Reallocate(newLength+1);

    copy(seq, seq + newLength, bases);
    length = newLength;
    bases[length] = '\0';
}

void Sequence::SetContent(const Kmer &kmer, int kmerLength)
{
    Clear();
    if (kmerLength+1 > capacity || kmerLength+1 < (capacity>>1))
        Reallocate(kmerLength+1);

    for (int i = 0; i < kmerLength; ++i)
        bases[i] = kmer.GetBase(i);
    length = kmerLength;
    bases[length] = '\0';
}

void Sequence::SetContent(CompactSequence *seq)
{
    Clear();
    if (seq->length+1 > capacity || seq->length+1 < (capacity>>1))
        Reallocate(seq->length+1);

    for (int i = 0; i < seq->length; ++i)
        bases[i] = seq->GetNucleotide(i);
    length = seq->length;
    bases[length] = '\0';
}

void Sequence::SetContent(Read *read)
{
    Clear();
    if (read->length+1 > capacity || read->length+1 < (capacity>>1))
        Reallocate(read->length+1);

    for (int i = 0; i < read->length; ++i)
        bases[i] = read->GetNucleotide(i);
    length = read->length;
    bases[length] = '\0';
}

bool Sequence::IsChar()
{
    for (int i = 0; i < length; ++i)
    {
        if (bases[i] != 'A' && bases[i] != 'C'
            && bases[i] != 'G' && bases[i] != 'T')
            return false;
    }

    return true;
}

bool Sequence::IsCodon()
{
    for (int i = 0; i < length; ++i)
    {
        if (bases[i] < 0 || bases[i] > 3)
            return false;
    }

    return true;
}



void Sequence::Encode()
{
    for (int i = 0; i < length; ++i)
    {
        bases[i] = code[(int)bases[i]];
//         switch (bases[i])
//         {
//             case 'A':
//                 bases[i] = 0;
//                 break;
//
//             case 'C':
//                 bases[i] = 1;
//                 break;
//
//             case 'G':
//                 bases[i] = 2;
//                 break;
//
//             case 'T':
//                 bases[i] = 3;
//                 break;
//
// #ifdef DEBUG
//             default:
//                 LogError("Sequence::Encode() encode failed\n");
//                 exit(1);
//                 break;
// #endif
//         }
    }
}

void Sequence::Decode()
{
    for (int i = 0; i < length; ++i)
    {
        bases[i] = code[(int)bases[i]];
//         switch(bases[i])
//         {
//             case 0:
//                 bases[i] = 'A';
//                 break;
//
//             case 1:
//                 bases[i] = 'C';
//                 break;
//
//             case 2:
//                 bases[i] = 'G';
//                 break;
//
//             case 3:
//                 bases[i] = 'T';
//                 break;
//
// #ifdef DEBUG
//             default:
//                 LogError("Sequence::Decode() decode failed\n");
//                 exit(1);
//                 break;
// #endif
//         }
    }
}

// void Sequence::Trim(int t)
// {
// #ifdef DEBUG
//     if (t > length)
//     {
//         LogError("Sequence::Trim() t is larger than length\n");
//         exit(1);
//     }
// #endif
//     length -= t;
// }

void Sequence::ReverseComplement()
{
    reverse(bases, bases + length);
    for (int i = 0; i < length; ++i)
    {
        switch(bases[i])
        {
            case 'A':
                bases[i] = 'T';
                break;

            case 'C':
                bases[i] = 'G';
                break;

            case 'G':
                bases[i] = 'C';
                break;

            case 'T':
                bases[i] = 'A';
                break;

            case 0:
            case 1:
            case 2:
            case 3:
                bases[i] = 3 - bases[i];
                break;

#ifdef DEBUG
            default:
                LogError("Sequence::ReverseComplement() failed\n");
                exit(1);
                break;
#endif
        }
    }
}

void Sequence::GetSubSequence(Sequence *seq, int offset, int subLength)
{
#ifdef DEBUG
    if (offset >= length || offset + subLength > length)
    {
        LogError("Sequence::GetSubSequence() out of range\n");
        exit(1);
    }
#endif
    seq->SetContent(bases + offset, subLength);
}

Kmer Sequence::GetKmer(int offset, int kmerLength)
{
#ifdef DEBUG
    if (offset + kmerLength > length)
    {
        LogError("GetKmer out of range\n");
        exit(1);
    }

    if (!IsCodon())
    {
        LogError("GetKmer from non-codon format\n");
        exit(1);
    }
#endif

    Kmer kmer = {{0}};
    for (int i = 0; i < kmerLength; ++i)
    {
        kmer.SetBase(i, bases[offset+i]);
    }
    return kmer;
}
