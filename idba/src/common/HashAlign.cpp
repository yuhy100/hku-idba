#include "globals.h"
#include "Contig.h"
#include "Utils.h"
#include "HashAlign.h"

#include <omp.h>

#include <iostream>
#include <algorithm>
#include <cstdio>

using namespace std;


void HashAlign::Initialize(vector<Contig> &contigs)
{
    int64 total = 0;
    for (int i = 0; i < (int)contigs.size(); ++i)
        total += contigs[i].Size();
    positions.reserve(total);

//#pragma omp parallel for
//    for (int i = 0; i < (int)contigs.size(); ++i)
//        InsertSequence((Sequence &)contigs[i]);

    //positions.resize(NumNodes());

//#pragma omp parallel for
    for (int i = 0; i < (int)contigs.size(); ++i)
    {
        InsertSequence(&contigs[i], i);
    }
}

void HashAlign::Initialize(vector<ContigNode> &contigs)
{
    int64 total = 0;
    for (int i = 0; i < (int)contigs.size(); ++i)
        total += contigs[i].GetSize();
    positions.reserve(total);

//#pragma omp parallel for
    for (int i = 0; i < (int)contigs.size(); ++i)
        InsertSequence(&contigs[i].GetContig(), i);
}

void HashAlign::InsertSequence(const Sequence &seq)
{
    if (seq.Size() < kmerLength)
        return;

    Kmer kmer;
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(seq[i]);
    for (int i = kmerLength-1; i < seq.Size(); ++i)
    {
        kmer.AddRight(seq[i]);
        KmerNode *p = InsertKmer(kmer);
        p->Data() = NumNodes() - 1;
    }
}

void HashAlign::InsertSequence(const Sequence *seq, int id)
{
    if (seq->Size() < kmerLength)
        return;

    Kmer kmer;
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight((*seq)[i]);
    for (int i = kmerLength-1; i < seq->Size(); ++i)
    {
        kmer.AddRight((*seq)[i]);

//        if (kmer.IsPalindrome())
//            continue;

        Kmer key = kmer;
        Kmer revComp = kmer;
        revComp.ReverseComplement();

        Position pos;
        pos.id = id;
        pos.offset = i - kmerLength + 1;
        pos.length = seq->Size();
        pos.isReverse = false;

        if (revComp < kmer)
        {
            key = revComp;
            pos.ReverseComplement();
        }

        KmerNode *p = InsertKmer(key);
        int index = positions.size();
        p->Data() = index;
        //positions[p->Data()] = pos;
        positions.push_back(pos);
        positions[index] = pos;
    }
}

void HashAlign::BuildTable(vector<Contig> &contigs)
{
}

void HashAlign::BuildPositions(vector<Contig> &contigs)
{
}

bool HashAlign::AlignKmer(const Kmer &kmer, Position &pos)
{
    KmerNode *p = GetNode(kmer);
    if (p == NULL)
        return false;

    pos = positions[p->Data()];
    if (kmer != p->GetKmer())
        pos.ReverseComplement();
    
    return true;
}

void HashAlign::AlignSequence(const Sequence *seq, vector<Alignment> &alignments)
{
    alignments.resize(0);

    if (seq->Size() < kmerLength)
        return ;
    {
        vector<Alignment> condidates;

        Kmer kmer;
        for (int i = 0; i < kmerLength-1; ++i)
            kmer.AddRight((*seq)[i]);
        for (int i = kmerLength-1; i < seq->Size(); ++i)
        {
            kmer.AddRight((*seq)[i]);

            Position pos;
            bool aligned = AlignKmer(kmer, pos);

            if (!aligned)
                continue;

            {
                Kmer revComp = kmer;
                revComp.ReverseComplement();

                bool flag = false;
                for (unsigned k = 0 ; k < condidates.size(); ++k)
                {
                    Alignment &align = condidates[k];
                    if (align.contigId == pos.id && align.isReverse == pos.isReverse 
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
                    align.readLength = seq->Size();
                    align.readOffset = i - kmerLength + 1;
                    align.contigId = pos.id;
                    align.contigOffset = pos.offset;
                    align.contigLength = pos.length;
                    align.length = kmerLength;
                    align.isReverse = pos.isReverse;
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
}

