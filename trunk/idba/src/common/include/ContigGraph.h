#ifndef __CONTIG_GRAPH_H_

#define __CONTIG_GRAPH_H_

#include "globals.h"
#include "Kmer.h"
#include "Contig.h"
#include "HashGraph.h"
#include "ContigNode.h"
#include "ContigBranchGroup.h"

#include <cstdio>
#include <vector>

class ContigGraph
{
public:
    ContigGraph()
    { hash_graph = new HashGraph(); }

    ContigGraph(std::vector<Contig> &contigs)
    { hash_graph = new HashGraph(); Initialize(contigs); }

    ~ContigGraph()
    { delete hash_graph; }

    void Clear() { hash_graph->Clear(); }

    void Initialize(std::vector<Contig> &contigs)
    { SetContigs(contigs); BuildVertices(); }
    void Initialize(std::vector<Contig> &contigs, std::vector<Kmer> &branches)
    { SetContigs(contigs), AddBranches(branches); BuildVertices(); }

    void BuildAllEdges()
    { hash_graph->AddAllEdges(); Refresh(); }
    bool BuildEdgesFromSequence(const Sequence &seq)
    { return hash_graph->AddEdgesFromSequence(seq); }

    void Refresh();
    bool Check();

    double AverageCoverage();

    int RemoveDeadEnd(int minLength);
    int RemoveLowCoverageContigs(double c);
    int RemoveBubble();
    void MergeContigs();

    int Assemble(std::vector<Contig> &contigs);
    int Assemble(std::vector<Contig> &contigs, std::vector<Kmer> &branches);

    ContigNodeAdapter GetNeighbor(const ContigNodeAdapter &node, int x);
    void GetNeighbors(const ContigNodeAdapter &node, std::vector<ContigNodeAdapter> &neighbors);

    std::vector<ContigNode> &GetContigNodes()
    { return nodes; }

protected:
    std::vector<ContigNode> nodes;

private:
    ContigGraph(const ContigGraph&);
    const ContigGraph &operator =(const ContigGraph &);

    void BuildVertices();
    void SetContigs(std::vector<Contig> &contigs);
    void AddBranches(std::vector<Kmer> &branches);

    //void RemoveContigNode(ContigNode &node);
    bool GetNextNodeAdapter(ContigNodeAdapter &current, ContigNodeAdapter &next);

    bool IsLoop(ContigPath &path, ContigNodeAdapter &next)
    {
        return path.GetBeginNodeAdapter().GetNode() == next.GetNode()
            || path.GetEndNodeAdapter().GetNode() == next.GetNode();
    }

    HashGraph *hash_graph;
};

#endif
