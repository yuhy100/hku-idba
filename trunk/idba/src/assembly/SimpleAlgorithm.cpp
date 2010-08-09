#include "globals.h"
#include "Log.h"
#include "Sequence.h"
#include "Utils.h"
#include "HashGraph.h"
#include "ContigGraph.h"
#include "HashAlign.h"
#include "MapGraph.h"
#include "SimpleAlgorithm.h"
#include "Reader.h"
#include "Writer.h"

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

SimpleAlgorithm::SimpleAlgorithm()
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


int64 SimpleAlgorithm::GenerateKmers()
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

void SimpleAlgorithm::AddbackInternalKmers()
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

void SimpleAlgorithm::RemoveLowCoverageContigs()
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

void SimpleAlgorithm::Iterate()
{
    ContigGraph contigGraph(contigs);

    LogMessage("start iteration from mink = %d to maxk = %d\n", kmerLength, maxk);
    while (kmerLength <= maxk)
    {
        contigGraph.Initialize(contigs);
        //contigGraph.BuildVertices();

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

        int minLength = kmerLength * 2;
        if (minLength < 75)
            minLength = 75;

        int deadend = 0;
        if (kmerLength < maxk)
        {
            deadend = contigGraph.RemoveDeadEnd(minLength);

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
            int bubbles = contigGraph.RemoveBubble();
            LogMessage("remove %d bubbles\n", bubbles);
        }

        contigGraph.MergeContigs();

        //int branches = contigGraph.BuildContigsFromEdges();
        //LogMessage("k = %d, remove %d deadend, remain %d branches\n", kmerLength, deadend, branches);

        ++kmerLength;
    }

    //contigGraph.GetContigs(contigs);
    contigGraph.Clear();
}

void SimpleAlgorithm::IDBA(FILE *freadFile, vector<Contig> &contigs)
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
    cout << cover << endl;

    hashGraph.Assemble(contigs);
    LogMessage("total kmer %d, Initial contigs %d\n", hashGraph.NumNodes(), contigs.size());
    hashGraph.Clear();

//    int tmp = 0;
//    for (int i = 0; i < contigs.size(); ++i)
//    {
//        if (contigs[i].Coverage() < cover)
//            ++tmp;
//    }
//    cout << tmp << endl;

    ContigGraph contigGraph(contigs);
    {
        contigGraph.Initialize(contigs);
        //contigGraph.BuildVertices();

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

        int minLength = kmerLength * 2;
        if (minLength < 75)
            minLength = 75;

//        {
//            int low = contigGraph.RemoveLowCoverageContigs(cover);
//            LogMessage("remove %d low coverage contigs\n", low);
//        }

        int deadend = 0;
        deadend = contigGraph.RemoveDeadEnd(minLength);
        LogMessage("remove %d deadend\n", deadend);
        contigGraph.MergeContigs();

        while (deadend > 0)
        {
            deadend = contigGraph.RemoveDeadEnd(minLength);
            LogMessage("remove %d deadend\n", deadend);
            contigGraph.MergeContigs();
        }

        int low = contigGraph.RemoveLowCoverageContigs(cover);
        LogMessage("remove %d low coverage contigs\n", low);
        contigGraph.MergeContigs();

        deadend = contigGraph.RemoveDeadEnd(minLength);
        LogMessage("remove %d deadend\n", deadend);
        contigGraph.MergeContigs();
        while (deadend > 0)
        {
            deadend = contigGraph.RemoveDeadEnd(minLength);
            LogMessage("remove %d deadend\n", deadend);
            contigGraph.MergeContigs();
        }

        contigGraph.MergeContigs();
    }

    //contigGraph.GetContigs(contigs);
    contigGraph.Clear();

    return;

    Iterate();

    if (kmerLength == maxk)
    {
        hashGraph.Clear();
#pragma omp parallel for
        for (int i = 0; i < (int)contigs.size(); ++i)
        {
            hashGraph.InsertSequence(contigs[i]);
        }

#pragma omp parallel for
        for (int i= 0; i < (int64)numReads; ++i)
        {
            Sequence seq;
            seq = reads[i];
            hashGraph.AddEdgesFromSequence(seq);
        }

#pragma omp parallel for
        for (int i = 0; i < (int)longReads.size(); ++i)
            hashGraph.AddEdgesFromSequence(longReads[i]);
        hashGraph.Refresh();

        LogMessage("k = %d, remain %lld kmers, %lld edges\n", kmerLength, hashGraph.NumNodes(), hashGraph.NumEdges());

        //hashGraph.RemoveDeadEnd(kmerLength*2);
        if (kmerLength == maxk)
        {
            int bubble = hashGraph.RemoveBubble();

            LogMessage("remove %d bubbles\n", bubble);
        }

        hashGraph.Assemble(contigs);
        hashGraph.Clear();
    }
    kmerLength = maxk;

    vector<Read>().swap(reads);
}


void SimpleAlgorithm::Run()
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
        writer.WriteFormat(contigs[i], "contig%d", i);//, i);
        contigs[i].Encode();
    }
}
