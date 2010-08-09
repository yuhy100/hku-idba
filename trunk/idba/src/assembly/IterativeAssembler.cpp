#include "globals.h"
#include "IterativeAssembler.h"
#include "Kmer.h"
#include "Contig.h"
#include "Log.h"

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <cmath>

using namespace std;

void IterativeAssembler::Assemble(vector<Contig> &contigs)
{
    FastAReader contig_reader(option.graphfile);

    Contig contig;
    string comment;
    while (contig_reader.Read(contig, comment))
    {
        int index = comment.size()-1;
        while (isdigit(comment[index]))
            --index;
        contig.sum_coverage = atoi(comment.c_str() + index + 1);
        contig.Encode();
        tmp_contigs.push_back(contig);
    }

    Iterate();

    contigs = tmp_contigs;
}

void IterativeAssembler::Iterate()
{
    vector<Kmer> branches;

    kmerLength = option.mink;
    LogMessage("start iteration from mink = %d to maxk = %d\n", kmerLength, option.maxk);
    while (kmerLength <= option.maxk)
    {
        contig_graph.Initialize(tmp_contigs, branches);

        RewindShortReads();
        vector<Read> &reads = GetShortReads();
#pragma omp parallel for
        for (int64 i = 0; i < (int64)reads.size(); ++i)
        {
            if (reads[i].IsActive() || reads[i].Size() < kmerLength)
            {
                Sequence seq;
                seq = reads[i];
                if (!contig_graph.BuildEdgesFromSequence(seq))
                    reads[i].Inactivate();
            }
        }

        vector<Sequence> &long_reads = GetLongReads();
#pragma omp parallel for
        for (int i = 0; i < (int)long_reads.size(); ++i)
            contig_graph.BuildEdgesFromSequence(long_reads[i]);
        contig_graph.Refresh();

        int minLength = kmerLength * 2;
        if (minLength < 75)
            minLength = 75;

        int deadend = 0;
        if (kmerLength < option.maxk)
        {
            if (kmerLength == option.mink)
            {
                option.cover = contig_graph.AverageCoverage();
                option.cover = sqrt(option.cover);
                if (option.cover < 2)
                    option.cover = 2;

                cout << option.cover << endl;
                //option.cover = 2;

                for (int i = 3; i >= 0; --i)
                {
                    int64 trimmed = contig_graph.RemoveDeadEnd(minLength >> i);
                    LogMessage("remove %d deadend\n", trimmed);
                    contig_graph.MergeContigs();
                }

                int64 low_coverage = contig_graph.RemoveLowCoverageContigs(option.cover);
                LogMessage("remove %d low coverage contigs\n", low_coverage);

                contig_graph.MergeContigs();
            }

            deadend = contig_graph.RemoveDeadEnd(minLength);
        }

        if (kmerLength == option.maxk)
        {
            int bubbles = contig_graph.RemoveBubble();
            LogMessage("remove %d bubbles\n", bubbles);
        }

        contig_graph.MergeContigs();
        int num = contig_graph.Assemble(tmp_contigs, branches);
        LogMessage("k = %d, remove %d deadend, remain %d branches\n", kmerLength, deadend, num);

        ++kmerLength;
    }

    contig_graph.Clear();
}
