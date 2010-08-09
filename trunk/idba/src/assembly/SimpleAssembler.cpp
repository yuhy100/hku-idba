#include "globals.h"
#include "Contig.h"
#include "SimpleAssembler.h"
#include "Log.h"

#include <vector>

using namespace std;

void SimpleAssembler::Assemble(vector<Contig> &contigs)
{
    kmerLength = option.mink;

    BuildKmerNodes();
    AddInternalKmerNodes();

    hash_graph->RefreshEdges();
    hash_graph->Assemble(contigs);
}

void SimpleAssembler::BuildKmerNodes()
{
    vector<vector<Kmer32> > solid_kmers(1 << option.prefix_length);

    int64 totalKmer = 0;
    int64 maximum = 0;

    for (int prefix = 0; prefix < (1 << option.prefix_length); ++prefix)
    {
        RewindShortReads();
        while (true)
        {
            vector<Read> &reads = GetShortReads();
            if (reads.size() == 0)
                break;


#pragma omp parallel for
            for (int i = 0; i < (int64)reads.size(); ++i)
            {
                if (!reads[i].IsActive() || reads[i].Size() < kmerLength)
                    continue;

                Sequence seq;
                seq = reads[i];
                hash_graph->InsertSequence(seq, prefix, option.mask);
            }
        }

        vector<Sequence> &long_reads = GetLongReads();
#pragma omp parallel for
        for (int i = 0; i < (int)long_reads.size(); ++i)
            hash_graph->InsertSequence(long_reads[i], prefix, option.mask);

        cout << "table size " << hash_graph->NumNodes() << " " << hash_graph->table_size << " " << hash_graph->backup.size()  << endl;
        hash_graph->RefreshVertices(option.min_count);

        hash_graph->Save(solid_kmers[prefix]);
        hash_graph->Clear();

        cout << "solid kmer " << solid_kmers[prefix].size() << endl;

        totalKmer = hash_graph->NumNodes();
        if (hash_graph->NumNodes() > maximum)
            maximum = hash_graph->NumNodes();
    }

    for (unsigned i = 0; i < solid_kmers.size(); ++i)
        hash_graph->Load(solid_kmers[i]);

    LogMessage("total kmer %d\n", hash_graph->NumNodes());
}

void SimpleAssembler::AddInternalKmerNodes()
{
    RewindShortReads();
    while (true)
    {
        vector<Read> &reads = GetShortReads();
        if (reads.size() == 0)
            break;

#pragma omp parallel for
        for (int i = 0; i < (int64)reads.size(); ++i)
        {
            if (!reads[i].IsActive() || reads[i].Size() < kmerLength)
                continue;

            Sequence seq;
            seq = reads[i];
            hash_graph->AddEdgesAndInternalKmers(seq);
        }
    }

    LogMessage("after addback total kmer %d\n", hash_graph->NumNodes());
}
