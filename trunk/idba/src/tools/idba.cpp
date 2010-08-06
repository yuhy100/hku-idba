#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "utils.h"
#include "hashGraph.h"
#include "compactSequence.h"
#include "read.h"

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

#include <pthread.h>

using namespace std;

const unsigned MaxThread = 4;

// static unsigned bits[] =
// {
//     0, 1, 1, 2, 1, 2, 2, 3,
//     1, 2, 2, 3, 2, 3, 3, 4,
// };

pthread_t threads[MaxThread];

unsigned minCount = 2;
int minLength = 0;
int maxK = 50;
int minK = kmerLength;
unsigned minContig = 100;
unsigned cover = 0;

uint64 numReads = 0;
int step = 1;
Sequence seq;
unsigned index2 = 0;
unsigned remove2 = 0;
unsigned round = 0;
unsigned oldMinCount = minCount;
int readLength = 0;

char line[MaxLine];
char comment[MaxLine];
unsigned char (*exist)[10];
Read *reads;
vector<Sequence> longReads;
vector<Sequence> contigs;
vector<Kmer> branches;
HashGraph *hashGraph = new HashGraph(MaxHashTable);

void *RefreshVerticesThread(void *p)
{
    uint64 i = (uint64)p;
    unsigned gap = hashGraph->tableSize / MaxThread;
    unsigned from = gap * i;
    unsigned to = gap * (i+1);
    if (i + 1 == MaxThread)
        to = hashGraph->tableSize;

    hashGraph->RefreshVertices(from, to, minCount);

    return p;
}

void *RefreshEdgesThread(void *p)
{
    uint64 i = (uint64)p;

    unsigned gap = hashGraph->tableSize / MaxThread;
    unsigned from = gap * i;
    unsigned to = gap * (i+1);
    if (i + 1 == MaxThread)
        to = hashGraph->tableSize;

    hashGraph->RefreshEdges(from, to);

    return p;
}

void *AddEdgesFromReadsThread(void *p)
{
    uint64 i = (uint64)p;
    unsigned gap = numReads / MaxThread;
    unsigned from = gap * i;
    unsigned to = gap * (i+1);
    if (i + 1 == MaxThread)
        to = numReads;

    Sequence seq;

    for (unsigned i = from; i < to; ++i)
    {
        if (reads[i].length >= kmerLength)
        {
            seq.SetContent(&reads[i]);
            int count = hashGraph->AddEdgesFromSequence(&seq);//, exist[i]);
            ++index2;

            if (count == 0)
            {
                reads[i].length = -reads[i].length;
                ++remove2;
            }
        }
    }

    return p;
}

void *RemoveDeadEndThread(void *p)
{
    uint64 i = (uint64)p;
    unsigned gap = hashGraph->tableSize / MaxThread;
    unsigned from = gap * i;
    unsigned to = gap * (i+1);
    if (i + 1 == MaxThread)
        to = hashGraph->tableSize;

    hashGraph->RemoveDeadEnd(minLength, from, to, HashGraphFlagScaning * i);

    return p;
}

void Refresh()
{
    for (unsigned i = 0; i < MaxThread; ++i)
        pthread_create(&threads[i], NULL, RefreshVerticesThread, (void *)i);

    for (unsigned i = 0; i < MaxThread; ++i)
        pthread_join(threads[i], NULL);

    for (unsigned i = 0; i < MaxThread; ++i)
        pthread_create(&threads[i], NULL, RefreshEdgesThread, (void *)i);

    for (unsigned i = 0; i < MaxThread; ++i)
        pthread_join(threads[i], NULL);
}

void AddEdgesFromReads()
{
    for (unsigned i = 0; i < MaxThread; ++i)
        pthread_create(&threads[i], NULL, AddEdgesFromReadsThread, (void *)i);

    for (unsigned i = 0; i < MaxThread; ++i)
        pthread_join(threads[i], NULL);

    for (unsigned i = 0; i < longReads.size(); ++i)
        hashGraph->AddEdgesFromSequence(&longReads[i]);
}

void RemoveDeadEnd()
{
    for (unsigned i = 0; i < MaxThread; ++i)
        pthread_create(&threads[i], NULL, RemoveDeadEndThread, (void *)i);

    for (unsigned i = 0; i < MaxThread; ++i)
        pthread_join(threads[i], NULL);
}

void InputReads(FILE *freadFile)
{
    numReads = CountSequences(freadFile);
    reads = new Read[numReads];
    exist = new unsigned char[numReads][10];
    fill_n(&exist[0][0], numReads * 10, (unsigned char)255);
    index2 = 0;
    while (ReadFasta(freadFile, &seq, comment))
    {
        if (seq.length >= (int)kmerLength && seq.IsChar())
        {
            seq.Encode();

            if (seq.length > (int)readLength)
                readLength = seq.length;

            if (seq.length <= 80)
                reads[index2++].SetContent(&seq);
            else
                longReads.push_back(seq);

            if (index2%2000000 == 0)
            {
                LogMessage("total reads %d\n", index2);
            }
        }
    }
    numReads = index2;
}

int main(int argc, char *argv[])
{
    kmerLength = 25;

    AddParameter("kmer", &kmerLength, INTEGER);
    AddParameter("minCount", &minCount, INTEGER);
//     AddParameter("minLength", &minLength, INTEGER);
//     AddParameter("step", &step, INTEGER);
    AddParameter("maxK", &maxK, INTEGER);
    AddParameter("minContig", &minContig, INTEGER);
    AddParameter("cover", &cover, INTEGER);

    ProcessParameters(argc, argv);

    if (minLength == 0)
        minLength = kmerLength * 2;

    if (minLength < readLength)
        minLength = readLength;

    minK = kmerLength;
    oldMinCount = minCount;

    if (argc < 4
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: idba read-file table-file contig-file\n");
        fprintf(stderr, "       [--kmer k] [--maxK k2] [--minCount m]\n");
        fprintf(stderr, "       [--minContig l] [--cover c]\n");
        exit(1);
    }

    FILE *freadFile = OpenFile(argv[1], "rb");
    FILE *ftableFile = OpenFile(argv[2], "rb");

    index2 = 0;

    ReadHashGraph(ftableFile, hashGraph, minCount);
    fclose(ftableFile);

    InputReads(freadFile);
    fclose(freadFile);

    if (readLength-1 < maxK)
        maxK = readLength - 1;

    remove2 = 0;
    round = 0;
    int oldN50 = 0;
    while (true)
    {
        LogMessage("Kmer %d maxK %d\n", kmerLength, maxK);

        if (round > 0)
        {
            hashGraph->Clear();

//             double start = clock();
            index2 = 0;
            for (unsigned i = 0; i < branches.size(); ++i)
            {
                hashGraph->InsertKmer(branches[i]);
            }
//             double end = clock();
//             fprintf(stderr, "Gen vertices from reads %.4f\n", (end - start)/CLOCKS_PER_SEC);

//             start = clock();
            hashGraph->ClearGraph();
            for (unsigned i = 0; i < contigs.size(); ++i)
            {
                if (contigs[i].length >= (int)kmerLength)
                {
                    contigs[i].Encode();
                    hashGraph->AddEdgesFromContig(&contigs[i]);
                }
            }

            index2 = 0;
            AddEdgesFromReads();
            for (unsigned i = 0; i < hashGraph->tableSize; ++i)
            {
                for (HashNode *node = hashGraph->table[i]; node; node = node->next)
                {
                    node->in |= node->inDepth;
                    node->out |= node->outDepth;
                    node->inDepth = node->outDepth = 0;
                }
            }

//             LogMessage("remove %d\n", remove2);
//             end = clock();
//             fprintf(stderr, "Gen edges from reads %.4f\n", (end - start)/CLOCKS_PER_SEC);
        }

        round++;

//         double start = clock();
        Refresh();
//         double end = clock();
//         fprintf(stderr, "refresh %.4f\n", (end - start)/CLOCKS_PER_SEC);
//         LogMessage("total kmer %d reads %d\n", hashGraph->numNodes, numReads);

        LogMessage("total reads %d/%d\ttotal kmer %d\t edges %d\n",
            index2, numReads, hashGraph->numNodes,
            hashGraph->numEdges);

//         start = clock();
        unsigned deadEnd = 0;
        if (round == 1)
        {
            RemoveDeadEnd();
            Refresh();

            if (round == 1)
            {
                int lowContigs = 0;
                if (cover == 0)
                    lowContigs = hashGraph->RemoveLowCoverageContigs(hashGraph->AverageCoverage()/5);
                else
                    lowContigs = hashGraph->RemoveLowCoverageContigs(cover);

                LogMessage("remove %d low coverage contigs\n", lowContigs);

                Refresh();
            }
        }

//         if (round > 1)
        if (kmerLength == maxK)
        {
            for (unsigned i = 0; i < hashGraph->tableSize; ++i)
            {
                for (HashNode *node = hashGraph->table[i]; node != NULL; node = node->next)
                {
                    node->count = 1;
                }
            }

            for (unsigned k = 0; k < numReads; ++k)
            {
                int l = reads[k].length;
                if (l < 0)
                    reads[k].length = -l;
//                 reads[k].length = readLength;

                if (reads[k].length >= (int)kmerLength)
                {
                    Sequence seq;
                    seq.SetContent(&reads[k]);

                    Kmer kmer = {{0}};
                    for (int i = 0; i < kmerLength-1; ++i)
                        kmer.AddRight(seq.bases[i]);
                    for (int i = kmerLength-1; i < seq.length; ++i)
                    {
                        kmer.AddRight(seq.bases[i]);
                        HashNode *node = hashGraph->GetByKmer(kmer);
                        if (node != NULL)
                        {
                            ++node->count;
                        }
                    }
                }

                reads[k].length = l;
            }

            RemoveDeadEnd();
            Refresh();

            int lowContigs = 0;
            if (cover == 0)
                lowContigs = hashGraph->RemoveLowCoverageContigs(hashGraph->AverageCoverage()/8);
            else
                lowContigs = hashGraph->RemoveLowCoverageContigs(cover);

            LogMessage("remove %d low coverage contigs\n", lowContigs);
            Refresh();
        }

//         end = clock();
//         fprintf(stderr, "check dead end %.4f\n", (end - start)/CLOCKS_PER_SEC);

//         LogMessage("Remove %d dead ends total %d \n", deadEnd, hashGraph->numNodes);
        unsigned remain = hashGraph->numNodes - deadEnd;
        LogMessage("remain %d kmers\n", remain);

        unsigned total = 0;
        if (kmerLength >= maxK)
        {
            unsigned bubble = hashGraph->RemoveBubble();
            LogMessage("Remove %d bubbles\n", bubble);
//             RemoveDeadEnd();
            Refresh();
            total = hashGraph->Assemble(contigs, branches);
        }
        else
        {
            if (kmerLength + step >= maxK)
                step = maxK - kmerLength;

//             start = clock();
            total = hashGraph->Assemble(contigs, branches);
//             end = clock();
//             fprintf(stderr, "assemble %.4f\n", (end - start)/CLOCKS_PER_SEC);
        }

        LogMessage("get contigs %d total nodes %d branches %llu\n",
            total, HashNode::totalNodes, branches.size());

        sprintf(line, "%s-%d", argv[3], kmerLength);
        FILE *fcontigFile = OpenFile(line, "wb");

        index2 = 0;
        int contigs100 = 0;
        for (unsigned i = 0; i < total; ++i)
        {
            if (contigs[i].length >= (int)minContig && contigs[i].length >= 2*(int)kmerLength)
                ++contigs100;

            {
                fprintf(fcontigFile, ">contig%d_length_%d\n",
                        index2++, contigs[i].length);
                WriteFasta(fcontigFile, &contigs[i]);
            }
        }
//         fprintf(stderr, "contigs100 %d\n", contigs100);
        fclose(fcontigFile);

        int *length = new int[contigs.size()];
        for (unsigned i = 0; i < total; ++i)
            length[i] = contigs[i].length;
        sort(length, length + total);
        reverse(length, length + total);

        int n50 = 0;
        long long sum = 0;
        for (int i = 0; i < contigs100; ++i)
        {
            sum += length[i];
            if (sum >= 0.5*remain && n50 == 0)
            {
                n50 = length[i];
            }
        }

        if (contigs100 == 0)
            contigs100 = 1;

        if (kmerLength >= maxK)
        {
            printf("read: %s \nkmer: %d minCount: %d maxK: %d \n",
                   argv[1], kmerLength, oldMinCount, maxK);
            printf("total: %d N50: %d ave: %lld max: %d nodes: %llu\n\n",
                   contigs100, n50, sum/contigs100, length[0], hashGraph->numNodes);
        }
        else
        {
            //fprintf(stderr, "kmer: %d\n", kmerLength);
            fprintf(stderr, "total: %d N50: %d ave: %lld max: %d\n\n",
                    contigs100, n50, sum/contigs100, length[0]);
        }

        delete [] length;

        if (kmerLength >= maxK)
            break;

        oldN50 = n50;
        if (kmerLength + step >= maxK)
            kmerLength = maxK;
        else
            kmerLength += step;

        minLength = kmerLength * 2;
        if (minLength < readLength)
            minLength = readLength;

        minCount = 1;
    }


    FILE *fcontigFile = OpenFile(argv[3], "wb");
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if (contigs[i].length >= (int)minContig && contigs[i].length >= 2*(int)kmerLength)
        {
            fprintf(fcontigFile, ">contig%d_length_%d\n",
                    index2++, contigs[i].length);
            WriteFasta(fcontigFile, &contigs[i]);
        }
    }
    fclose(fcontigFile);

    return 0;
}
