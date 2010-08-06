#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "utils.h"
#include "read.h"
#include "hashTable.h"

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;

struct Position
{
    int id;
    int offset;

    void ReverseComplement();
};

struct Alignment
{
    int readLength;
    int readOffset;
    int contigId;
    int contigOffset;
    int contigLength;
    int length;
    int isReverse;
};

int totalContigs;
int alignedReads;
char line[MaxLine];
char comment[MaxLine];
HashTable *hashTable = new HashTable(MaxHashTable, true);
Sequence *contigs;
Position *positions;
Position posBuffer[MaxLine];

void Position::ReverseComplement()
{
    offset = contigs[id>>1].length - (offset + kmerLength);
    id ^= 1;
}

void PreAddSequence(const Sequence *seq)
{
    Kmer kmer = {{0}};
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(seq->bases[i]);
    for (int i = kmerLength-1; i < seq->length; ++i)
    {
        kmer.AddRight(seq->bases[i]);
        hashTable->InsertKmer(kmer);

        Kmer revCmp = kmer;
        revCmp.ReverseComplement();
        if (revCmp == kmer)
            hashTable->InsertKmer(revCmp);
    }
}

void AllocatePositions()
{
    unsigned index = 0;
    for (unsigned i = 0; i < hashTable->tableSize; ++i)
    {
        for (HashNode *p = hashTable->table[i]; p != NULL; p = p->next)
        {
            p->Data() = index;
            index += p->count;
            p->count = 0;
        }

//         if (i%100000 == 0)
//             cout << i << " " << index << endl;
    }
    positions = new Position[index];
}

void AddSequence(const Sequence *seq, int id)
{
    Kmer kmer = {{0}};
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(seq->bases[i]);
    for (int i = kmerLength-1; i < seq->length; ++i)
    {
        kmer.AddRight(seq->bases[i]);
        HashNode *p = hashTable->InsertKmer(kmer);
        int index = p->Data() + p->count - 1;

        positions[index].id = id;
        positions[index].offset = i - kmerLength + 1;
        if (kmer != p->kmer)
            positions[index].ReverseComplement();

        Kmer revCmp = kmer;
        revCmp.ReverseComplement();
        if (revCmp == kmer)
        {
            HashNode *p = hashTable->InsertKmer(revCmp);
            int index = p->Data() + p->count - 1;

            positions[index].id = id;
            positions[index].offset = i - kmerLength + 1;
            positions[index].ReverseComplement();
        }
    }
}

int GetPosBuffer(const Kmer &kmer)
{
    HashNode *p = hashTable->GetByKmer(kmer);
    int count = p->count;
    int index = p->Data();
    for (int i = 0; i < count; ++i)
    {
        posBuffer[i] = positions[index + i];
        if (kmer != p->kmer)
            posBuffer[i].ReverseComplement();
//         posBuffer[k+1] = posBuffer[k];
//         posBuffer[k+1].ReverseComplement();
//
//         if (kmer != p->kmer)
//             swap(posBuffer[k], posBuffer[k+1]);
    }

    return count;
}

int main(int argc, char *argv[])
{
    AddParameter("kmer", &kmerLength, INTEGER);

    ProcessParameters(argc, argv);

    if (argc < 4
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: alignedReads read-file contig-file align-file\n");
        fprintf(stderr, "       [--kmer k]\n");
        exit(1);
    }

    FILE *freadFile = OpenFile(argv[1], "rb");
    FILE *fcontigs = OpenFile(argv[2], "rb");
    FILE *falign = OpenFile(argv[3], "wb");

    totalContigs = CountSequences(fcontigs);
    contigs = new Sequence[totalContigs];

    int index = 0;
    Sequence seq;
    while (ReadFasta(fcontigs, &seq, comment))
    {
        seq.Encode();
        PreAddSequence(&seq);
        contigs[index++] = seq;
    }

//     cout << hashTable->numNodes << endl;

    AllocatePositions();

    for (int i = 0; i < totalContigs; ++i)
        AddSequence(&contigs[i], i << 1);

//      cout << hashTable->numNodes << endl;

    index = 0;
    while (ReadFasta(freadFile, &seq, comment))
    {
        if (seq.IsChar())
        {
            seq.Encode();

            int hits = 0;
            vector<Alignment> alignments;

            {
                vector<Alignment> condidates;

                Kmer kmer = {{0}};
                for (int i = 0; i < kmerLength-1; ++i)
                    kmer.AddRight(seq.bases[i]);
                for (int i = kmerLength-1; i < seq.length; ++i)
                {
                    kmer.AddRight(seq.bases[i]);
                    HashNode *p = hashTable->GetByKmer(kmer);

                    if (p == NULL)
                        continue;

                    ++hits;
                    int numPos = GetPosBuffer(kmer);
                    for (int j = 0; j < numPos; ++j)
                    {
                        Position pos = posBuffer[j];

                        bool flag = false;
                        for (unsigned k = 0 ; k < condidates.size(); ++k)
                        {
                            Alignment &align = condidates[k];
                            if (align.contigId == pos.id
                                && align.contigOffset + align.length == pos.offset + kmerLength - 1)
                            {
                                align.length += 1;
                                flag = true;
                                break;
                            }
                        }

                        if (!flag)
                        {
                            Alignment align;
                            align.readLength = seq.length;
                            align.readOffset = i - kmerLength + 1;
                            align.contigId = pos.id;
                            align.contigOffset = pos.offset;
                            align.contigLength = contigs[pos.id>>1].length;
                            align.length = kmerLength;
                            condidates.push_back(align);
                        }
                    }

                    int l = 0;
                    for (unsigned k = 0; k < condidates.size(); ++k)
                    {
                        Alignment &align = condidates[k];
                        if (i < align.readOffset + align.length)
                            condidates[l++] = align;
                        else
                            alignments.push_back(align);
                    }
                    condidates.resize(l);
                }

                for (unsigned k = 0; k < condidates.size(); ++k)
                    alignments.push_back(condidates[k]);

            }

            fprintf(falign, "read %d alignments %d hits %d\n",
                    index, (int)alignments.size(), hits);
            seq.Decode();

            for (unsigned i = 0; i < alignments.size(); ++i)
            {
                Alignment &align = alignments[i];

                fprintf(falign, "readLength %d readoffset %d contigId %d contigOffset %d contigLength %d length %d reverse %d\n",
                        align.readLength, align.readOffset,
                        align.contigId, align.contigOffset, align.contigLength,
                        align.length, align.isReverse);
            }

            if (alignments.size() > 0)
            {
                ++alignedReads;
            }
        }
        else
        {
            fprintf(falign, "read %d alignments %d hits %d\n",
                    index, 0, 0);
        }

        index++;
    }


    fprintf(falign, ">kmer: %d total contigs: %d total kmers: %llu\n",
            kmerLength, totalContigs, hashTable->numNodes);
    fprintf(falign, ">alignedReads %d/%d\n", alignedReads, index);

    return 0;
}