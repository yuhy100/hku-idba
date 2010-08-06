#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "utils.h"
#include "read.h"
#include "hashTable.h"
#include "contigGraph.h"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;

struct Alignment
{
    int readLength;
    int readOffset;
    int contigId;
    int contigOffset;
    int contigLength;
    int length;
    int isReverse;

    void ReverseComplement()
    {
        readOffset = readLength - (readOffset + length);
        contigOffset = contigLength - (contigOffset + length);
        contigId ^= 1;
        isReverse = !isReverse;
    }
};


char line[MaxLine];
char comment[MaxLine];
char buf[MaxLine];
int numReads;
int minPairs = 10;
int hits[MaxDistance] = {0};
int offset = 10;
ContigGraph contigGraph;
int estimateDistance;
int estimateRange;
int minContig = 100;

bool ReadAlignment(FILE *fp, vector<Alignment> &alignments)
{
    int index;
    int n;

    if (fgets(line, MaxLine, fp) == NULL)
        return false;

    if (line[0] == '>')
        return false;

    sscanf(line, "%s %d %s %d", buf, &index, buf, &n);
    Alignment align;
    for (int i = 0; i < n; ++i)
    {
        fgets(line, MaxLine, fp);
        sscanf(line, "%s %d %s %d %s %d %s %d %s %d %s %d %s %d",
               buf, &align.readLength,
               buf, &align.readOffset,
               buf, &align.contigId,
               buf, &align.contigOffset,
               buf, &align.contigLength,
               buf, &align.length,
               buf, &align.isReverse);
        if (align.length >= kmerLength && align.length >= 0.6*align.readLength)
            alignments.push_back(align);
    }

    return true;
}

void EstimateDistance(FILE *fp)
{
    double sumDistance = 0;
    int numDistance = 0;
    int maxDistance = 0;
    while (true)
    {
        vector<Alignment> u;
        vector<Alignment> v;

        if (!ReadAlignment(fp, u) || !ReadAlignment(fp, v))
            break;

        for (unsigned i = 0; i < u.size(); ++i)
        {
            for (unsigned j = 0; j < v.size(); ++j)
            {
                Alignment a = u[i];
                Alignment b = v[j];

                b.ReverseComplement();

                if (a.contigId == b.contigId)
                {
                    int from = a.contigOffset - a.readOffset;
                    int to = b.contigOffset + (b.readLength - b.readOffset);
                    if (to-from < 0)
                        continue;

                    ++numDistance;
                    if (to - from < (int)MaxDistance)
                    {
                        ++hits[to-from];
                        if (to - from > maxDistance)
                        {
                            maxDistance = to - from;
                        }
                    }
                }
            }
        }
    }

    int discard = numDistance/200;
    int from = 0;
    int to = maxDistance;
    int sum = 0;
    while (sum + hits[from] < discard)
    {
        sum += hits[from++];
    }

    sum = 0;
    while (sum + hits[to] < discard)
    {
        sum += hits[to--];
    }

    sumDistance = 0;
    int realNum = 0;
    for (int i = from; i <= to; ++i)
    {
        sumDistance += i * hits[i];
        realNum += hits[i];
    }

    estimateDistance = int(round(sumDistance/realNum));

    int region = numDistance * 90 / 100;
    sum = hits[estimateDistance];
    offset = 1;
    while (offset < estimateDistance && sum  < region)
    {
        sum += hits[estimateDistance - offset] + hits[estimateDistance + offset];
        ++offset;
    }

    if (offset < estimateDistance/20)
        offset = estimateDistance/20;

    estimateRange = offset;
    printf("estimate distance %d +- %d max %d\n",
           estimateDistance, estimateRange, maxDistance);

    fseek(fp, 0, SEEK_SET);
}

int main(int argc, char *argv[])
{
    AddParameter("kmer", &kmerLength, INTEGER);
    AddParameter("minPairs", &minPairs, INTEGER);
    AddParameter("minContig", &minContig, INTEGER);
    ProcessParameters(argc, argv);

    if (argc < 4
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: mergedContigs contig-file align-file mate-contig-file\n");
        fprintf(stderr, "       [--kmer k] [--minPairs n] [--minContig l]\n");
        exit(1);
    }

    FILE *fcontigs = OpenFile(argv[1], "rb");
    FILE *falign = OpenFile(argv[2], "rb");
    FILE *fout = OpenFile(argv[3], "wb");

    contigGraph.ReadFromFile(fcontigs);
    fclose(fcontigs);

    EstimateDistance(falign);

    while (true)
    {
        vector<Alignment> u;
        vector<Alignment> v;

        if (!ReadAlignment(falign, u) || !ReadAlignment(falign, v))
            break;

        if (u.size() > 2 && v.size() > 2)
            continue;

        for (unsigned i = 0; i < u.size(); ++i)
        {
            for (unsigned j = 0; j < v.size(); ++j)
            {
                Alignment a = u[i];
                Alignment b = v[j];

                if (a.contigId > b.contigId)
                    swap(a, b);

                b.ReverseComplement();

                if (a.contigId != b.contigId)
                {
                    int x = a.contigId;
                    int y = b.contigId;
                    int d = estimateDistance - (a.contigLength - (a.contigOffset - a.readOffset))
                            - (b.contigOffset + (b.readLength - b.readOffset));

                    if (d < -2*kmerLength || d > estimateDistance)
                        continue;

                    contigGraph.AddConnection(x, y, d);

                    swap(x, y);
                    x ^= 1;
                    y ^= 1;

                    contigGraph.AddConnection(x, y, d);
                }
            }
        }
    }

    contigGraph.SetRange(estimateRange);
    contigGraph.FilterConnections(minPairs);
    contigGraph.BuildInitialPaths();
    contigGraph.FindPossibleConnections();

    vector<Sequence> mergedContigs;
    contigGraph.MergeContigs(mergedContigs);

    int index = 0;
    for (unsigned i = 0; i < mergedContigs.size(); ++i)
    {
        if (mergedContigs[i].length >= minContig)
        {
            fprintf(fout, ">contig%d_length%d\n", index++, mergedContigs[i].length);
            WriteFasta(fout, &mergedContigs[i]);
        }
    }

    fclose(falign);
    fclose(fout);

    return 0;
}
