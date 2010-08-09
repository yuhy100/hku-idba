#include "globals.h"
#include "GappedScaffoldAssembler.h"

#include <cstdio>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

void GappedScaffoldAssembler::Assemble(vector<Contig> &contigs)
{
    kmerLength = option.maxk;
    scaffold_graph.SetMinPairs(option.min_pairs);
    FastAReader contig_reader(option.scaffile);

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

    cout << "aligned" << endl;
    scaffold_graph.ComputeDistance();
    scaffold_graph.ScaffoldWithGap(contigs);
}
