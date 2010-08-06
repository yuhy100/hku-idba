#include "globals.h"
#include "sequence.h"
#include "kmer.h"
#include "hashNode.h"
#include "hashGraph.h"
#include "log.h"

#include <cstdio>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <queue>

using namespace std;

typedef HashNode *HashNodePointer;

static const uint64 UnitOne = 100000000ULL;

static unsigned bits[] =
{
    0, 1, 1, 2, 1, 2, 2, 3,
    1, 2, 2, 3, 2, 3, 3, 4,
};

static unsigned bitToIndex[] =
{
    0, 0, 1, 0, 2, 0, 0, 0, 3,
};

HashGraph::HashGraph(uint64 tableSize)
{
    table = NULL;
    this->tableSize = 0;
    numNodes = 0;
    numEdges = 0;

    Reallocate(tableSize);
}

HashGraph::~HashGraph()
{
    Clear();
    delete [] table;
}

void HashGraph::Reallocate(uint64 newTableSize)
{
    HashNode **newTable = new HashNodePointer[newTableSize];

    if (newTable == NULL)
    {
        LogError("HashGraph::Reallocate() not enough memory\n");
        exit(1);
    }

    for (uint64 i = 0; i < newTableSize; ++i)
        newTable[i] = NULL;

    for (uint64 i = 0; i < tableSize; ++i)
    {
        HashNode *node = table[i];
        while (node != NULL)
        {
            HashNode *next = node->next;
            uint64 index = node->kmer.Hash() % newTableSize;
            node->next = newTable[index];
            newTable[index] = node;
            node = next;
        }
    }

    table = newTable;
    tableSize = newTableSize;
}

void HashGraph::Clear()
{
    numNodes = 0;
    numEdges = 0;
    for (unsigned i = 0; i < tableSize; ++i)
    {
        HashNode *node = table[i];
        while (node != NULL)
        {
            HashNode *p = node;
            node = node->next;
            HashNode::FreeNode(p);
        }
        table[i] = NULL;
    }
}

void HashGraph::ClearGraph()
{
    for (unsigned i = 0; i < tableSize; ++i)
    {
        for (HashNode *node = table[i]; node != NULL; node = node->next)
        {
            node->status = 0;
            node->inDepth = 0;
            node->outDepth = 0;
            node->count = 1;
            node->in = 0;
            node->out = 0;
        }
    }

    numEdges = 0;
}

HashNode *HashGraph::InsertKmer(const Kmer &kmer, int count)
{
    if (numNodes > tableSize)
        Reallocate(tableSize * 2);

    Kmer key = kmer;
    Kmer revComp = kmer;
    revComp.ReverseComplement();
    if (revComp < kmer)
        key = revComp;

    uint64 index = key.Hash() % tableSize;
    HashNode *p = table[index];

    while (p != NULL && p->kmer != key)
        p = p->next;

    if (p == NULL)
    {
        p = HashNode::NewNode();
        p->SetContent(key);
        p->next = table[index];
        table[index] = p;
        ++numNodes;
    }

    if (count + p->count < 65535)
        p->count += count;
    else
        p->count = 65535;

    return p;
}

HashNode *HashGraph::GetByKmer(const Kmer &kmer)
{
    Kmer key = kmer;
    Kmer revComp = kmer;
    revComp.ReverseComplement();
    if (revComp < kmer)
        key = revComp;

    uint64 index = key.Hash() % tableSize;
    HashNode *p = table[index];

    while (p != NULL && p->kmer != key)
        p = p->next;

    return p;
}

void HashGraph::InsertSequence(const Sequence *seq)
{
    Kmer kmer = {{0}};
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(seq->bases[i]);
    for (int i = kmerLength-1; i < seq->length; ++i)
    {
        kmer.AddRight(seq->bases[i]);
        InsertKmer(kmer);
    }
}

bool HashGraph::IsValid(Sequence *seq)
{
    Kmer kmer = {{0}};
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(seq->bases[i]);
    for (int i = kmerLength-1; i < seq->length; ++i)
    {
        kmer.AddRight(seq->bases[i]);
        if (GetByKmer(kmer) == NULL)
            return false;
    }

    return true;
}

int HashGraph::AddEdgesFromSequence(const Sequence *seq)
{
    int branch = 0;
    int count = 0;
    Kmer kmer = {{0}};
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(seq->bases[i]);

    Kmer x = {{0}};
    HashNode *last = NULL;
    for (int i = kmerLength-1; i < seq->length; ++i)
    {
        kmer.AddRight(seq->bases[i]);

        HashNode *node = GetByKmer(kmer);
        if (node != NULL)
        {
            node->count++;
            ++count;

            if (last != NULL)
            {
                HashNodeAdapter adp(last, x);
                HashNodeAdapter adp2(node, kmer);

                if ((adp.OutEdges() & (1 << seq->bases[i])) == 0
                      || (adp2.InEdges() & (1 << (3 - seq->bases[i-kmerLength]))) == 0)
                    ++branch;

                adp.OutDepth() |= 1 << seq->bases[i];
                adp2.InDepth() |= 1 << (3 - seq->bases[i-kmerLength]);
            }
        }

        x = kmer;
        last = node;
    }

    if (branch == 0)
        return 0;

    return count;
}

int HashGraph::AddEdgesFromSequence(const Sequence *seq, unsigned char *exist)
{
    unsigned char flags[10] = {0};

    int branch = 0;
    int count = 0;
    Kmer kmer = {{0}};
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(seq->bases[i]);

    Kmer x = {{0}};
    HashNode *last = NULL;
//     unsigned char *prevOut = NULL;
//     unsigned char *prevOutDepth = NULL;
    for (int i = kmerLength-1; i < seq->length; ++i)
    {
        kmer.AddRight(seq->bases[i]);

        HashNode *node = NULL;
        if (exist[i>>3] & (1 << (i&7)))
        {
            node = GetByKmer(kmer);
            if (node != NULL)
            {
                node->count++;
                ++count;

//                 unsigned char *currIn = &node->in;
//                 unsigned char *currOut = &node->out;
//                 unsigned char *currInDepth = &node->inDepth;
//                 unsigned char *currOutDepth = &node->outDepth;
//                 if (node->kmer != kmer)
//                 {
//                     swap(currIn, currOut);
//                     swap(currInDepth, currOutDepth);
//                 }
//
//                 if (prevOut != NULL)
//                 {
//                     flags[i>>3] |= 1 << (i & 7);
//                     if ((*prevOut & (1 << seq->bases[i])) == 0
//                           || (*currIn & (1 << (3 - seq->bases[i-kmerLength]))) == 0)
//                         ++branch;
//
//                     *prevOutDepth |= 1 << seq->bases[i];
//                     *currInDepth |= 1 << (3 - seq->bases[i-kmerLength]);
//                 }

                if (last != NULL)
                {
                    flags[i>>3] |= 1 << (i & 7);
                    HashNodeAdapter adp(last, x);
                    HashNodeAdapter adp2(node, kmer);

                    if ((adp.OutEdges() & (1 << seq->bases[i])) == 0
                         || (adp2.InEdges() & (1 << (3 - seq->bases[i-kmerLength]))) == 0)
                        ++branch;

                    adp.OutDepth() |= 1 << seq->bases[i];
                    adp2.InDepth() |= 1 << (3 - seq->bases[i-kmerLength]);
                }

//                 prevOut = currOut;
//                 prevOutDepth = currOutDepth;
            }
//             else
//                 prevOut = NULL;
        }
//         else
//             prevOut = NULL;

        x = kmer;
        last = node;
    }

    if (branch == 0)
        return 0;

    copy(flags, flags + 10, exist);

    return count;
}

int HashGraph::AddEdgesFromContig(const Sequence *contig)
{
    int count = 0;
    Kmer kmer = {{0}};
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(contig->bases[i]);

    Kmer x = {{0}};
    HashNode *last = NULL;
    for (int i = kmerLength-1; i < contig->length; ++i)
    {
        kmer.AddRight(contig->bases[i]);

        HashNode *node = InsertKmer(kmer);
        if (node != NULL)
        {
            node->count++;
            ++count;

            if (last != NULL)
            {
                HashNodeAdapter adp(last, x);
                HashNodeAdapter adp2(node, kmer);

                adp.OutEdges() |= 1 << contig->bases[i];
                adp2.InEdges() |= 1 << (3 - contig->bases[i-kmerLength]);
            }
        }

        x = kmer;
        last = node;
    }

    return count;
}

uint64 HashGraph::RemoveDeadEnd(unsigned minLength, unsigned from, unsigned to, unsigned scan)
{
    if (from == 0 && to == 0)
        to = tableSize;

    unsigned count = 0;
    queue<Kmer> qu;
    for (unsigned i = from; i < to; ++i)
    {
        for (HashNode *node = table[i]; node != NULL; node = node->next)
        {
            if (node->GetStatus(HashGraphFlagAlive))
                continue;

            Kmer kmer = node->kmer;

            if (CheckDeadEnd(kmer, minLength, scan))
            {
                node->SetStatus(HashGraphFlagDeadend);
                count += 1;
                continue;
            }
            kmer.ReverseComplement();
            qu.push(kmer);
            while (!qu.empty())
            {
                Kmer k = qu.front();
                qu.pop();

                HashNode *p = GetByKmer(k);
                HashNodeAdapter adapter(p, k);

                if (adapter.InDepth() < minLength)
                {
                    adapter.InDepth() = minLength;
                    unsigned out = adapter.OutEdges();
                    for (int x = 0; x < 4; ++x)
                    {
                        if (out & (1 << x))
                        {
                            qu.push(k);
                            qu.back().AddRight(x);
                        }
                    }
                }
            }

            if (CheckDeadEnd(kmer, minLength, scan))
            {
                node->SetStatus(HashGraphFlagDeadend);
                count += 1;
                continue;
            }
            kmer.ReverseComplement();
            qu.push(kmer);
            while (!qu.empty())
            {
                Kmer k = qu.front();
                qu.pop();

                HashNode *p = GetByKmer(k);
                HashNodeAdapter adapter(p, k);

                if (adapter.InDepth() < minLength)
                {
                    adapter.InDepth() = minLength;
                    unsigned out = adapter.OutEdges();
                    for (int x = 0; x < 4; ++x)
                    {
                        if (out & (1 << x))
                        {
                            qu.push(k);
                            qu.back().AddRight(x);
                        }
                    }
                }
            }

            node->SetStatus(HashGraphFlagAlive);
        }
    }

    LogMessage("remove deadend %d\n", count);

    unsigned recue = 0;
    for (unsigned i = from; i < to; ++i)
    {
        for (HashNode *node = table[i]; node != NULL; node = node->next)
        {
            if (node->GetStatus(HashGraphFlagDeadend))
                continue;

            Kmer kmer = node->kmer;
            for (int strand = 0; strand < 2; ++strand)
            {
                HashNode *current = GetByKmer(kmer);
                HashNodeAdapter adapter(current, kmer);
                unsigned out = adapter.OutEdges();

                Kmer x = kmer;
                while (true)
                {
                    if (bits[out] == 1)
                    {
                        x.AddRight(bitToIndex[out]);
                        HashNode *next = GetByKmer(x);
                        HashNodeAdapter adapter2(next, x);

                        if (next != current
                            && next->GetStatus(HashGraphFlagDeadend)
                            && bits[adapter2.OutEdges()] == 1 && bits[adapter2.InEdges()] == 1)
                        {
                            ++recue;
                            next->SetStatus(HashGraphFlagAlive);
                            next->ResetStatus(HashGraphFlagDeadend);
                            current = next;
                            out = adapter2.OutEdges();
                        }
                        else
                            break;
                    }
                    else
                        break;
                }

                kmer.ReverseComplement();
            }
        }
    }

    LogMessage("rescue %d\n", recue);

    return count - recue;
}

int HashGraph::CheckDeadEnd(const Kmer &kmer, unsigned remain, unsigned scan)
{
    if (remain == 0)
    {
        return false;
    }
    else
    {
        HashNode *node = GetByKmer(kmer);
        HashNodeAdapter adapter(node, kmer);

        if (node->GetStatus(scan))
            return true;

        if (node->GetStatus(HashGraphFlagAlive))
            return false;

        if (adapter.OutDepth() >= remain)
            return false;

        unsigned edges = adapter.OutEdges();
        bool flag = true;
        node->SetStatus(scan);

        for (int i = 0; i < 4; ++i)
        {
            if (edges & (1 << i))
            {
                Kmer x = kmer;
                x.AddRight(i);
                if (!CheckDeadEnd(x, remain-1, scan))
                {
                    flag = false;
                    break;
                }
            }
        }

        node->ResetStatus(scan);

        return flag;
    }
}

void HashGraph::Refresh(unsigned minCount)
{
    numEdges = 0;
    for (unsigned i = 0; i < tableSize; ++i)
    {
        HashNode *node = table[i];
        HashNode *prev = NULL;
        while (node != NULL)
        {
            if (node->GetStatus(HashGraphFlagDeadend) || node->count < minCount)
            {
                --numNodes;
                if (prev == NULL)
                {
                    table[i] = node->next;
                    HashNode::FreeNode(node);
                    node = table[i];
                }
                else
                {
                    prev->next = node->next;
                    HashNode::FreeNode(node);
                    node = prev->next;
                }
            }
            else
            {
                HashNodeAdapter adapter(node);
                Kmer kmer = adapter.node->kmer;
                for (int strand = 0; strand < 2; ++strand)
                {
                    unsigned edges = adapter.OutEdges();

                    for (int x = 0; x < 4; ++x)
                    {
                        if (edges & (1 << x))
                        {
                            Kmer next = kmer;
                            next.AddRight(x);
                            HashNode *q = GetByKmer(next);
                            if (q == NULL || q->GetStatus(HashGraphFlagDeadend)
                               || q->count < minCount)
                            {
                                adapter.OutEdges() ^= (1 << x);
                            }
                            else
                            {
                                ++numEdges;
                            }
                        }
                    }

                    adapter.ReverseComplement();
                    kmer.ReverseComplement();
                }

                prev = node;
                node = prev->next;
            }
        }
    }

    numEdges >>= 1;
}

void HashGraph::RefreshVertices(unsigned from, unsigned to, unsigned minCount)
{
    if (from == 0 && to == 0)
    {
        to = tableSize;
    }

    for (unsigned i = from; i < to; ++i)
    {
        HashNode *node = table[i];
        HashNode *prev = NULL;
        while (node != NULL)
        {
            bool dead = false;
            if (node->GetStatus(HashGraphFlagDeadend) || node->count < minCount)
                dead = true;

            if (!dead)
            {
                HashNodeAdapter adp(node);
                Kmer tmp = adp.node->kmer;
                for (int strand = 0; strand < 2; ++strand)
                {
                    if (bits[adp.InEdges()] == 0 && bits[adp.OutEdges()] == 1)
                    {
                        tmp.AddRight(bitToIndex[adp.OutEdges()]);
                        HashNode *p = GetByKmer(tmp);
                        if (p != NULL)
                        {
                            HashNodeAdapter adp2(p, tmp);
                            if (bits[adp2.InEdges()] > 1)
                            {
                                dead = true;
                                break;
                            }
                        }
                    }

                    adp.ReverseComplement();
                    tmp.ReverseComplement();
                }

                node->ResetStatus(HashGraphFlagAlive);
                node->inDepth = node->outDepth = 0;
            }

            if (dead)
            {
                --numNodes;
                if (prev == NULL)
                {
                    table[i] = node->next;
                    HashNode::FreeNode(node);
                    node = table[i];
                }
                else
                {
                    prev->next = node->next;
                    HashNode::FreeNode(node);
                    node = prev->next;
                }
            }
            else
            {
                prev = node;
                node = prev->next;
            }
        }
    }
}

void HashGraph::RefreshEdges(unsigned from, unsigned to)
{
    numEdges = 0;
    if (from == 0 && to == 0)
    {
        to = tableSize;
    }

    for (unsigned i = from; i < to; ++i)
    {
        HashNode *node = table[i];
        HashNode *prev = NULL;
        while (node != NULL)
        {
            HashNodeAdapter adapter(node);
            Kmer kmer = adapter.node->kmer;
            for (int strand = 0; strand < 2; ++strand)
            {
                unsigned edges = adapter.OutEdges();

                for (int x = 0; x < 4; ++x)
                {
                    if (edges & (1 << x))
                    {
                        Kmer next = kmer;
                        next.AddRight(x);
                        HashNode *q = GetByKmer(next);
                        if (q == NULL || q->GetStatus(HashGraphFlagDeadend))
                        {
                            adapter.OutEdges() ^= (1 << x);
                        }
                        else
                        {
                            ++numEdges;
                        }
                    }
                }

                kmer.ReverseComplement();
                adapter.ReverseComplement();
            }

            Kmer revComp = node->kmer;
            revComp.ReverseComplement();
            if (revComp == node->kmer)
                node->in = node->out = node->in | node->out;

            prev = node;
            node = prev->next;
        }
    }

    numEdges >>= 1;
}

bool HashGraph::ErrorCorrect(Sequence *seq, int index, Kmer kmer, unsigned error)
{
    if (index == seq->length)
        return true;
    else
    {
        HashNode *node = GetByKmer(kmer);

        if (node == NULL || node->GetStatus(HashGraphFlagDeadend))
            return false;

        HashNodeAdapter adp(node, kmer);
        unsigned out = adp.OutEdges();
        {
            unsigned c = seq->bases[index];

            if (out & (1 << c))
            {
                Kmer x = kmer;
                x.AddRight(c);
                if (ErrorCorrect(seq, index+1, x, error))
                    return true;
            }

            for (unsigned i = 0; i < 4; ++i)
            {
                if ((out & (1 << i)) && c != i && error > 0)
                {
                    seq->bases[index] = i;
                    Kmer x = kmer;
                    x.AddRight(i);
                    if (ErrorCorrect(seq, index+1, x, error - 1))
                        return true;
                }
            }

            seq->bases[index] = c;
        }

        return false;
    }
}

void HashGraph::ErrorCorrectWeight(Sequence *seq, int index, Kmer kmer, unsigned error, int sum)
{
    if (index == seq->length)
    {

        HashNode *node = GetByKmer(kmer);

        if (node == NULL || node->GetStatus(HashGraphFlagDeadend))
            return;

        sum += node->count;

        if (sum > correctedWeight)
        {
            if (correctedWeight >= 0)
            {
//                 cerr << correctedWeight << " " << sum << endl;
                Sequence seq1 = corrected;
                Sequence seq2 = *seq;

                seq1.Decode();
                seq2.Decode();

//                 cerr << seq1 << endl;
//                 cerr << seq2 << endl;
            }

            correctedWeight = sum;
            corrected = *seq;
        }
    }
    else
    {
        HashNode *node = GetByKmer(kmer);

        if (node == NULL || node->GetStatus(HashGraphFlagDeadend))
            return;

        HashNodeAdapter adp(node, kmer);
        unsigned out = adp.OutEdges();
        {
            unsigned c = seq->bases[index];

            if (out & (1 << c))
            {
                Kmer x = kmer;
                x.AddRight(c);

                ErrorCorrectWeight(seq, index+1, x, error, sum + node->count);
            }
            //else
            {
                for (unsigned i = 0; i < 4; ++i)
                {
                    if ((out & (1 << i)) && c != i && error > 0)
                    {
                        seq->bases[index] = i;
                        Kmer x = kmer;
                        x.AddRight(i);

                        ErrorCorrectWeight(seq, index+1, x, error - 1, sum + node->count);
                    }
                }
            }

            seq->bases[index] = c;
        }
    }
}

uint64 HashGraph::RemoveBubble()
{
//     LogMessage("nodes %d edges %d\n", numNodes, numEdges);
    unsigned cand = 0;
    int index = 0;
    unsigned bubble = 0;
    for (unsigned i = 0; i < tableSize; ++i)
    {
        for (HashNode *node = table[i]; node != NULL; node = node->next)
        {
            ++index;

            HashNodeAdapter adp(node);
            if (!((bits[adp.OutEdges()] == 2 && bits[adp.InEdges()] == 1)
                   || (bits[adp.OutEdges()] == 1 && bits[adp.InEdges()] == 2)))
                continue;

            ++cand;

            Kmer kmer = node->kmer;
            for (int strand = 0; strand < 2; ++strand)
            {
                if (bits[adp.OutEdges()] == 2)
                {
                    vector<Kmer> path1;
                    vector<Kmer> path2;

                    uint64 weight1 = UnitOne / node->count;
                    uint64 weight2 = UnitOne / node->count;

                    path1.push_back(kmer);
                    for (int j = 0; j <= kmerLength; ++j)
                    {
                        Kmer from = path1[j];
                        HashNode *p = GetByKmer(from);
                        HashNodeAdapter adp2(p, from);
                        unsigned edges = adp2.OutEdges();

                        Kmer to = from;
                        if (j == 0)
                        {
                            to.AddRight(bitToIndex[edges & ~(edges - 1)]);
                        }
                        else if (bits[edges] == 1)
                        {
                            to.AddRight(bitToIndex[edges]);
                        }
                        else
                            break;

                        HashNode *q = GetByKmer(to);
                        HashNodeAdapter adp3(q, to);
                        if (q->GetStatus(HashGraphFlagDeadend))
                            break;
                        else if ((j == kmerLength && bits[adp3.InEdges()] == 2)
                                  || (j < kmerLength && bits[adp3.InEdges()] == 1))
                        {
                            weight1 += UnitOne / q->count;
                            path1.push_back(to);
                        }
                        else
                            break;
                    }

                    path2.push_back(kmer);
                    for (int j = 0; j <= kmerLength; ++j)
                    {
                        Kmer from = path2[j];
                        HashNode *p = GetByKmer(from);
                        HashNodeAdapter adp2(p, from);
                        unsigned edges = adp2.OutEdges();

                        Kmer to = from;
                        if (j == 0)
                        {
                            to.AddRight(bitToIndex[edges & (edges - 1)]);
                        }
                        else if (bits[edges] == 1)
                        {
                            to.AddRight(bitToIndex[edges]);
                        }
                        else
                            break;

                        HashNode *q = GetByKmer(to);
                        HashNodeAdapter adp3(q, to);
                        if (q->GetStatus(HashGraphFlagDeadend))
                            break;
                        else if ((j == kmerLength && bits[adp3.InEdges()] == 2)
                                  || (j < kmerLength && bits[adp3.InEdges()] == 1))
                        {
                            weight2 += UnitOne / q->count;
                            path2.push_back(to);
                        }
                        else
                            break;
                    }

                    HashNode *start = GetByKmer(path1[0]);
                    HashNodeAdapter adp2(start, path1[0]);
                    HashNode *end = GetByKmer(path1.back());
                    HashNodeAdapter adp3(end, path1.back());
                    if ((int)path1.size() == kmerLength + 2
                        && (int)path2.size() == kmerLength + 2
                        && path1.back() == path2.back()
                        && start != end
                        && bits[adp2.InEdges()] == 1
                        && bits[adp3.OutEdges()] == 1)
                    {
                        ++bubble;

                        if (weight1 > weight2)
                            path1.swap(path2);

                        adp2.OutEdges() &= (1 << path1[1].GetBase(kmerLength-1));
                        adp3.InEdges() &= (1 << (3 - path1[kmerLength].GetBase(0)));

                        if (bits[adp2.OutEdges()] != 1 || bits[adp3.InEdges()] != 1)
                        {
                            cerr << "bubble error" << endl;
                        }

                        for (int j = 1; j <= kmerLength; ++j)
                        {
                            HashNode *p = GetByKmer(path2[j]);
                            p->SetStatus(HashGraphFlagDeadend);
                        }

                    }
                }

                kmer.ReverseComplement();
                adp.ReverseComplement();
            }
        }
    }

//     LogMessage("candidate %d\n", cand);

    return bubble;
}

// uint64 HashGraph::RemoveLowCoverageContigs()
// {
//     long long sum = 0;
//     uint64 valid = 0;
//     for (uint64 i = 0; i < tableSize; ++i)
//     {
//         for (HashNode *node = table[i]; node; node = node->next)
//         {
//             sum += node->count;
//             ++valid;
//         }
//     }
//
//     double ave = 1.0 * sum / valid;
//     LogMessage("sum %lld ave %.4f\n", sum, ave);
//
//     return RemoveLowCoverageContigs(ave/5);
// }

double HashGraph::AverageCoverage()
{
    long long sum = 0;
    uint64 valid = 0;
    for (uint64 i = 0; i < tableSize; ++i)
    {
        for (HashNode *node = table[i]; node; node = node->next)
        {
            sum += node->count;
            ++valid;
        }
    }

    return 1.0 * sum / valid;
}

uint64 HashGraph::RemoveLowCoverageContigs(double c)
{
    unsigned total = 0;

    for (unsigned i = 0; i < tableSize; ++i)
    {
        for (HashNode *node = table[i]; node != NULL; node = node->next)
        {
            if (bits[node->in] > 1 || bits[node->out] > 1)
                node->SetStatus(HashGraphFlagUsed);
            Kmer kmer = node->kmer;
            kmer.ReverseComplement();

            if (kmer == node->kmer)
                node->SetStatus(HashGraphFlagUsed);
        }
    }

    for (unsigned i = 0; i < tableSize; ++i)
    {
        for (HashNode *node = table[i]; node != NULL; node = node->next)
        {
            if ((node->status & HashGraphFlagUsed)
                 || (node->status & HashGraphFlagDeadend))
            {
                continue;
            }

            vector<HashNode *> path;
            path.push_back(node);

            node->status |= HashGraphFlagUsed;

            Sequence contig;
            contig.SetContent(node->kmer, kmerLength);

            long long sum = node->count;
            int length = 1;

            for (int strand = 0; strand < 2; ++strand)
            {
                Kmer kmer =
                        contig.GetKmer(contig.length - kmerLength, kmerLength);

                while (true)
                {
                    HashNode *p = GetByKmer(kmer);
                    HashNodeAdapter adp(p, kmer);
                    unsigned edges = adp.OutEdges();


                    if (bits[edges] == 1)
                    {
                        int x = bitToIndex[edges];

                        Kmer next = kmer;
                        next.AddRight(x);
                        HashNode *q = GetByKmer(next);
                        HashNodeAdapter adp2(q, next);

                        if (bits[adp2.InEdges()] != 1
                            || (q->status & HashGraphFlagUsed)
                            || (q->status & HashGraphFlagDeadend))
                        {
                            break;
                        }

                        q->status |= HashGraphFlagUsed;
                        sum += q->count;
                        ++length;
                        path.push_back(q);

                        contig.AddNucleotide(x);
                        kmer = next;
                    }
                    else
                        break;
                }

                contig.ReverseComplement();
            }

            contig.Decode();


            if (length > 10 && 1.0*sum/length < c)
            {
                total++;

                for (unsigned i = 0; i < path.size(); ++i)
                    path[i]->SetStatus(HashGraphFlagDeadend);
            }
        }
    }

    for (unsigned i = 0; i < tableSize; ++i)
    {
        for (HashNode *node = table[i]; node != NULL; node = node->next)
        {
            node->ResetStatus(HashGraphFlagUsed);
        }
    }

    return total;
}

uint64 HashGraph::Assemble(Sequence *contigs)
{
    unsigned total = 0;

    for (unsigned i = 0; i < tableSize; ++i)
    {
        for (HashNode *node = table[i]; node != NULL; node = node->next)
        {
            if ((node->status & HashGraphFlagUsed)
                 || (node->status & HashGraphFlagDeadend))
            {
                continue;
            }

            node->status |= HashGraphFlagUsed;

            Sequence contig;
            contig.SetContent(node->kmer, kmerLength);

            for (int strand = 0; strand < 2; ++strand)
            {
                Kmer kmer =
                        contig.GetKmer(contig.length - kmerLength, kmerLength);

                while (true)
                {
                    HashNode *p = GetByKmer(kmer);
                    HashNodeAdapter adp(p, kmer);
                    unsigned edges = adp.OutEdges();


                    if (bits[edges] == 1)
                    {
                        int x = bitToIndex[edges];

                        Kmer next = kmer;
                        next.AddRight(x);
                        HashNode *q = GetByKmer(next);
                        HashNodeAdapter adp2(q, next);

                        if (bits[adp2.InEdges()] != 1
                            || (q->status & HashGraphFlagUsed)
                            || (q->status & HashGraphFlagDeadend))
                        {
                            break;
                        }

                        q->status |= HashGraphFlagUsed;

                        contig.AddNucleotide(x);
                        kmer = next;
                    }
                    else
                        break;
                }

                contig.ReverseComplement();
            }

            contig.Decode();
            contigs[total].SetContent(contig.bases, contig.length);
            total++;
        }
    }

    return total;
}

uint64 HashGraph::Assemble(vector<Sequence> &contigs, vector<Kmer> &branches)
{
    branches.resize(0);
    contigs.resize(0);
    unsigned total = 0;
    Sequence seq;
    Sequence contig;
    int tangle = 0;

    for (unsigned i = 0; i < tableSize; ++i)
    {
        for (HashNode *node = table[i]; node != NULL; node = node->next)
        {
            if ((node->status & HashGraphFlagUsed)
                 || (node->status & HashGraphFlagDeadend))
            {
                continue;
            }

            node->status |= HashGraphFlagUsed;

            contig.SetContent(node->kmer, kmerLength);
            bool end = false;

            for (int strand = 0; strand < 2; ++strand)
            {
                Kmer kmer =
                        contig.GetKmer(contig.length - kmerLength, kmerLength);
                unsigned edges = 0;

                while (true)
                {
                    HashNode *p = GetByKmer(kmer);
                    HashNodeAdapter adp(p, kmer);
                    edges = adp.OutEdges();

                    if (bits[edges] == 1)
                    {
                        int x = bitToIndex[edges];
                        Kmer next = kmer;
                        next.AddRight(x);
                        HashNode *q = GetByKmer(next);
                        HashNodeAdapter adp2(q, next);

                        if (q->status & HashGraphFlagDeadend)
                            break;
                        else if (bits[adp2.InEdges()] != 1)
                        {
                            break;
                        }
                        else if (q->status & HashGraphFlagUsed)
                        {
                            break;
                        }
                        else
                        {
                            q->status |= HashGraphFlagUsed;
                            contig.AddNucleotide(x);
                            kmer = next;
                        }
                    }
                    else
                    {
                        if (bits[edges] == 0)
                        {
                            end = true;
                        }

                        break;
                    }

                }

                for (int i = 0; i < 4; ++i)
                {
                    if (edges & (1 << i))
                    {
                        Kmer x = kmer;
                        x.SetBase(kmerLength, i);
                        branches.push_back(x);
                    }
                }

                contig.ReverseComplement();
            }

            if (end && contig.length < (int)kmerLength * 3)
                continue;

            if (end)
                ++tangle;

            contig.Decode();
            contigs.push_back(contig);
            total++;
        }
    }

//     cerr << branches.size() << endl;
    LogMessage("tangle %d total %d\n", tangle, total);
//     cerr << " tangle " << tangle << " total " << total << endl;

    return total;
}
