#include "globals.h"
#include "ScaffoldAssembler.h"
#include "Reader.h"

#include <cstdio>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

void ScaffoldAssembler::Assemble(vector<Contig> &contigs)
{
    kmerLength = option.maxk;
    scaffold_graph.SetMinPairs(option.min_pairs);
    FastAReader contig_reader(option.contigfile);

    Contig contig;
    string comment;
    while (contig_reader.Read(contig, comment))
    {
        int index = comment.size()-1;
        while (isdigit(comment[index]))
            --index;
        contig.sum_coverage = atoi(comment.c_str() + index + 1);
        contig.Encode();
        contigs.push_back(contig);
    }

    scaffold_graph.Initialize(contigs);

    cout << "scaffold graph initialize" << endl;
    aligner.Initialize(scaffold_graph.GetContigNodes());
    cout << "aligner initialize" << endl;

    Sequence seq1, seq2;

    while (true)
    {
        vector<Read> reads = GetShortReads();
        if (reads.size() == 0)
            break;

        for (int64 i = 0; i < (int64)reads.size(); i += 2)
        {
            if (!reads[i].IsActive() || !reads[i+1].IsActive())
                continue;

            vector<Alignment> alignments1;
            vector<Alignment> alignments2;
            seq1 = reads[i];
            seq2 = reads[i+1];

            aligner.AlignSequence(&seq1, alignments1);
            aligner.AlignSequence(&seq2, alignments2);

            scaffold_graph.AddPair(alignments1, alignments2);
        }
    }
    scaffold_graph.ComputeDistance();

    if (option.readfile2 != "")
    {
        FastAReader reader(option.readfile2);

        string comment1;
        string comment2;
        while (reader.Read(seq1, comment1) && reader.Read(seq2, comment2))
        {
            if (!seq1.IsChar() || !seq2.IsChar())
                continue;

            seq1.Encode();
            seq2.Encode();

            vector<Alignment> alignments1;
            vector<Alignment> alignments2;

            aligner.AlignSequence(&seq1, alignments1);
            aligner.AlignSequence(&seq2, alignments2);

            scaffold_graph.AddPair(alignments1, alignments2, 1);
        }

        cout << "read2 aligned" << endl;
        scaffold_graph.ComputeDistance(1);
        cout << "ok" << endl;
    }

    cout << "aligned" << endl;
    scaffold_graph.FindUniquePaths();
    scaffold_graph.Scaffold(contigs);
}
