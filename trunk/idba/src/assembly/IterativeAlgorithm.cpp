#include "globals.h"
#include "Log.h"
#include "Sequence.h"
#include "Utils.h"
#include "HashGraph.h"
#include "ContigGraph.h"
#include "HashAlign.h"
#include "MapGraph.h"
#include "IterativeAlgorithm.h"
#include "Reader.h"
#include "Writer.h"
#include "Kmer.h"

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

IterativeAlgorithm::IterativeAlgorithm()
{
    prefixLength = 3;
    minCount = 2;
    trim = 0;
    cover = 0;
    isScaffold = false;

    minPairs = 5;
    fill_n(hits, MaxDistance, 0);
    offset = 10;
    minContig = 100;

    prefix = "out";
    mink = 25;
    maxk = 50;
    minPairs = 5;
    minLength = 100;
}


int64 IterativeAlgorithm::GenerateKmers()
{
    int64 totalKmer = 0;
    int64 maximum = 0;

    for (int prefix = 0; prefix < (1 << prefixLength); ++prefix)
    {
#pragma omp parallel for
        for (int i = 0; i < (int64)numReads; ++i)
        {
            Sequence seq;
            seq = reads[i];
            hashGraph.InsertSequence(seq, prefix, mask);
        }

#pragma omp parallel for
        for (int i = 0; i < (int)longReads.size(); ++i)
        {
            hashGraph.InsertSequence(longReads[i], prefix, mask);
        }

        hashGraph.RefreshVertices(minCount);

        totalKmer = hashGraph.NumNodes();
        if (hashGraph.NumNodes() > maximum)
            maximum = hashGraph.NumNodes();

        //hashGraph.Save(caches[prefix]);
    }

    LogMessage("total kmer %d\n", hashGraph.NumNodes());

    return hashGraph.NumNodes();
}

void IterativeAlgorithm::AddbackInternalKmers()
{
#pragma omp parallel for
    for (int i = 0; i < (int64)numReads; ++i)
    {
        Sequence seq;
        //seq.SetContent(&reads[i]);
        seq = reads[i];
        hashGraph.AddInternalKmers(seq, minCount);
    }

    LogMessage("after addback total kmer %d\n", hashGraph.NumNodes());
}

void IterativeAlgorithm::RemoveLowCoverageContigs()
{
    if (cover == 0)
    {
        cover = hashGraph.MedianCoverage();
        cover = sqrt(cover);
        if (cover < 2)
            cover = 2;
        int lowCoverageContigs = hashGraph.RemoveLowCoverageContigs(cover);
        LogMessage("c = %.2f remove %d low coverage contigs\n", cover, lowCoverageContigs);
    }
    else
        hashGraph.RemoveLowCoverageContigs(cover);
    hashGraph.Refresh();
}

void IterativeAlgorithm::Iterate()
{
    vector<Kmer> branches;
    ContigGraph contigGraph(contigs);

    LogMessage("start iteration from mink = %d to maxk = %d\n", kmerLength, maxk);
    while (kmerLength <= maxk)
    {
        //contigGraph.GetContigs(contigs);

        //contigGraph.SetContigs(contigs);
        //if (kmerLength < maxk)
        //contigGraph.AddBranches(branches);
        contigGraph.Initialize(contigs, branches);

        //contigGraph.BuildVertices();
        //contigGraph.Refresh();
        //cout << "Build vertices ok" << endl;

#pragma omp parallel for
        for (int64 i = 0; i < (int64)numReads; ++i)
        {
            if (reads[i].IsActive())
            {
                Sequence seq;
                seq = reads[i];
                if (!contigGraph.BuildEdgesFromSequence(seq))
                    reads[i].Inactivate();
            }
        }

#pragma omp parallel for
        for (int i = 0; i < (int)longReads.size(); ++i)
            contigGraph.BuildEdgesFromSequence(longReads[i]);
        contigGraph.Refresh();
        //cout << "Build contig graph ok" << endl;

        int minLength = kmerLength * 2;
        if (minLength < 75)
            minLength = 75;

        int deadend = 0;
        if (kmerLength < maxk)
        {
            if (kmerLength == mink)
            {
                cout << cover << endl;
                cover = contigGraph.AverageCoverage();
                cover = sqrt(cover);
                cout << cover << endl;
                if (cover < 2)
                    cover = 2;

                for (int i = 3; i >= 0; --i)
                {
                    int64 trimmed = contigGraph.RemoveDeadEnd(minLength >> i);
                    LogMessage("remove %d deadend\n", trimmed);
                    contigGraph.MergeContigs();
                }

                int64 low_coverage = contigGraph.RemoveLowCoverageContigs(cover);
                LogMessage("remove %d low coverage contigs\n", low_coverage);

                contigGraph.MergeContigs();
            }

            deadend = contigGraph.RemoveDeadEnd(minLength);
            //cout << "remove deadend ok" << endl;

//            if (kmerLength == mink)
//            {
//                LogMessage("remove %d deadend\n", deadend);
////                int low = contigGraph.RemoveLowCoverageContigs(cover);
////                LogMessage("remove %d low coverage contigs\n", low);
//                while (deadend > 0)
//                {
//                    deadend = contigGraph.RemoveDeadEnd(minLength);
//                    LogMessage("remove %d deadend\n", deadend);
//                }
//            }
        }
        //contigGraph.Refresh();

        if (kmerLength == maxk)
        {
            //contigGraph.MergeContigs();
            int bubbles = contigGraph.RemoveBubble();
            LogMessage("remove %d bubbles\n", bubbles);
        }

        contigGraph.MergeContigs();
        //cout << "merge contigs ok" << endl;

        //int branches = contigGraph.BuildContigsFromEdges();
        //int branches = contigGraph.Assemble(contigs, branches);
        int num = contigGraph.Assemble(contigs, branches);
        LogMessage("k = %d, remove %d deadend, remain %d branches\n", kmerLength, deadend, num);

        ++kmerLength;

        //contigGraph.GetContigs(contigs);
    }

    contigGraph.Clear();
}

void IterativeAlgorithm::IDBA(FILE *freadFile, vector<Contig> &contigs)
{
    LogMessage("minK = %d, # of reads = %lld\n", kmerLength, numReads);

    GenerateKmers();
    AddbackInternalKmers();
    Sequence seq;

    hashGraph.RefreshEdges();
//    hashGraph.RemoveDeadEnd(kmerLength*2);
//    RemoveLowCoverageContigs();
    
    cover = hashGraph.MedianCoverage();
    cover = sqrt(cover);
    if (cover < 2)
        cover = 2;

    hashGraph.Assemble(contigs);
    LogMessage("total kmer %d, Initial contigs %d\n", hashGraph.NumNodes(), contigs.size());
    hashGraph.Clear();

    Iterate();
    kmerLength = maxk;

    vector<Read>().swap(reads);
}

void IterativeAlgorithm::EstimateDistance()
{
    double sumDistance = 0;
    int numDistance = 0;
    int maxDistance = 0;
    for (int i = 0; i < (int)alignResults.size(); i += 2)
    {
        vector<Alignment> &u = alignResults[i];
        vector<Alignment> &v = alignResults[i+1];

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
}

void IterativeAlgorithm::Align(const Sequence &seq, vector<Alignment> &alignments)
{
    alignments.resize(0);
    hashAlign.AlignSequence(&seq, alignments);

    for (unsigned i = 0; i < alignments.size(); ++i)
    {
        Alignment &align = alignments[i];
        align.contigId = (align.contigId << 1) + (!!align.isReverse);
    }

    int l = 0;
    for (unsigned i = 0; i < alignments.size(); ++i)
    {
        Alignment &align = alignments[i];
        if (align.length >= seq.Size()/2)
        {
            alignments[l++] = align;
        }
    }
    alignments.resize(l);
}

//void IterativeAlgorithm::Align(FILE *freadFile, vector<Contig> &contigs)
void IterativeAlgorithm::Align(Reader &reader, vector<Contig> &contigs)
{
    totalContigs = contigs.size();
    hashAlign.Initialize(contigs);

    reader.Rewind();
    //fseek(freadFile, 0, SEEK_SET);
    Sequence seq;
    int index = 0;
    alignResults.reserve(numReads);
    string comment;
    //while (ReadFasta(freadFile, &seq, comment))
    while (reader.Read(seq, comment))
    {
        seq.Trim(trim);
        seq.TrimError();
        if (seq.IsChar() && seq.Size() >= kmerLength)
        {
            seq.Encode();

            vector<Alignment> alignments;
            Align(seq, alignments);
//            hashAlign.AlignSequence(&seq, alignments);
//
//            seq.Decode();
//
//            for (unsigned i = 0; i < alignments.size(); ++i)
//            {
//                Alignment &align = alignments[i];
//                align.contigId = (align.contigId << 1) + (!!align.isReverse);
//            }
//
//            int l = 0;
//            for (unsigned i = 0; i < alignments.size(); ++i)
//            {
//                Alignment &align = alignments[i];
//                if (align.length >= seq.Size()/2)
//                {
//                    alignments[l++] = align;
//                }
//            }
//            alignments.resize(l);

            alignResults.push_back(alignments);

            if (alignments.size() > 0)
                ++alignedReads;
        }
        else
        {
            vector<Alignment> alignments;
            alignResults.push_back(alignments);
        }

        index++;
    }

    hashAlign.Clear();

    LogMessage("aligned %d reads\n", alignedReads);
}

void IterativeAlgorithm::MergeContigs(vector<Sequence> &mergedContigs)
{
    for (int i = 0; i < (int)alignResults.size(); i += 2)
    {
        vector<Alignment> &u = alignResults[i];
        vector<Alignment> &v = alignResults[i+1];

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

                    mapGraph.AddConnection(x, y, d);

                    swap(x, y);
                    x ^= 1;
                    y ^= 1;

                    mapGraph.AddConnection(x, y, d);
                }
            }
        }
    }

    mapGraph.SetRange(estimateRange);
    mapGraph.BuildInitialPaths();
    mapGraph.FilterConnections(minPairs);
    mapGraph.FindPossibleConnections();

    mapGraph.MergeContigs(mergedContigs);
}

//void IterativeAlgorithm::MergeContigs2(FILE *freadFile, std::vector<Contig> &contigs, std::vector<Sequence> &mergedContigs)
void IterativeAlgorithm::MergeContigs2(Reader &reader, std::vector<Contig> &contigs, std::vector<Sequence> &mergedContigs)
{
    totalContigs = contigs.size();
    hashAlign.Initialize(contigs);

//    while (true)
//    {
//        cout << "hello2" << endl;
//    }
//
    double sumDistance = 0;
    int numDistance = 0;
    int maxDistance = 0;

    reader.Rewind();
    //fseek(freadFile, 0, SEEK_SET);
    Sequence seq;
    Sequence seq2;
    string comment;
    string comment2;
    vector<int> x_cache;
    vector<int> y_cache;
    vector<int> d_cache;
    //while (ReadFasta(freadFile, &seq, comment) && ReadFasta(freadFile, &seq2, comment2))
    while (reader.Read(seq, comment) && reader.Read(seq2, comment2))// &seq, comment) && ReadFasta(freadFile, &seq2, comment2))
    {
        seq.Trim(trim);
        seq.TrimError();
        seq2.Trim(trim);
        seq2.TrimError();

        if (seq.IsChar() && seq.Size() >= kmerLength && seq2.IsChar() && seq2.Size() >= kmerLength)
        {
            seq.Encode();
            seq2.Encode();

            vector<Alignment> u;
            vector<Alignment> v;
            Align(seq, u);
            Align(seq2, v);

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
                        int d = - (a.contigLength - (a.contigOffset - a.readOffset))
                                - (b.contigOffset + (b.readLength - b.readOffset));

                        //if (d < -2*kmerLength || d > estimateDistance)
                        //    continue;

                        //mapGraph.AddConnection(x, y, d);

                        x_cache.push_back(x);
                        y_cache.push_back(y);
                        d_cache.push_back(d);
                            
                        swap(x, y);
                        x ^= 1;
                        y ^= 1;

                        x_cache.push_back(x);
                        y_cache.push_back(y);
                        d_cache.push_back(d);
                        //mapGraph.AddConnection(x, y, d);
                    }
                }
            }
        }
    }

    hashAlign.Clear();

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

    for (unsigned i = 0; i < x_cache.size(); ++i)
    {
        int x = x_cache[i];
        int y = y_cache[i];
        int d = estimateDistance + d_cache[i];
        if (d < -2*kmerLength || d > estimateDistance)
            continue;
        mapGraph.AddConnection(x, y, d);
    }

    mapGraph.SetDistance(estimateDistance);
    mapGraph.SetRange(estimateRange);
    mapGraph.BuildInitialPaths();
    mapGraph.FilterConnections(minPairs);
    mapGraph.FindPossibleConnections();

    mapGraph.MergeContigs(mergedContigs);
}

void IterativeAlgorithm::Run()
{
    minContig = minLength;
    kmerLength = mink;
    contigfile = prefix + "-contig.fa";
    scaffile = prefix + "-contig-mate.fa";

    mask = (1 << prefixLength) - 1;
    FILE *freadFile = OpenFile(readfile.c_str(), "rb");

    numReads = LoadReads(readfile);

    IDBA(freadFile, contigs);


    FastAWriter writer(contigfile);
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        contigs[i].Decode();
        writer.WriteFormat(contigs[i], "contig%d_%d", i, contigs[i].Size());
        contigs[i].Encode();
    }


    if (!isScaffold)
        return;

    kmerLength = maxk;
    mapGraph.Initialize(contigs);

    vector<Sequence> mergedContigs;
    FastAReader reader(readfile);
    MergeContigs2(reader, contigs, mergedContigs);

    FastAWriter writer2(scaffile);
    int index = 0;
    for (unsigned i = 0; i < mergedContigs.size(); ++i)
    {
        //if (mergedContigs[i].Size() >= minContig)
        {
            writer2.WriteFormat(mergedContigs[i], "contig%d_length%d", index++, mergedContigs[i].Size());
        }
    }

}
