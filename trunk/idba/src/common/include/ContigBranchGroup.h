#ifndef __CONTIG_BRANCH_GROUP_H_

#define __CONTIG_BRANCH_GROUP_H_

#include "globals.h"
#include "Contig.h"
#include "ContigNode.h"

#include <vector>
#include <algorithm>

class ContigGraph;

class ContigPath
{
public:
    ContigPath() { weight = 0; }

    void Clear() { weight = 0; nodes.resize(0); }

    ContigNodeAdapter &GetBeginNodeAdapter()
    { return nodes.front(); }
    ContigNodeAdapter &GetEndNodeAdapter()
    { return nodes.back(); }

    ContigNodeAdapter &operator [](int index) 
    { return nodes[index]; }

    uint64 Weight() { return weight; }
    int64 Size() { return nodes.size(); }
    int64 InternalSize();
    int64 InternalDistance();

    void Append(const ContigNodeAdapter &node)
    { nodes.push_back(node); weight += UnitOne / node.GetNode()->GetContig().Coverage(); }

    void Pop()
    { weight -= UnitOne / nodes.back().GetNode()->GetContig().Coverage(); nodes.pop_back(); }

    void ReverseComplement()
    {
        reverse(nodes.begin(), nodes.end());
        for (unsigned i = 0; i < nodes.size(); ++i)
            nodes[i].ReverseComplement();
    }

    bool IsSimplePath();
    void Inactivate();
    void Activate();

    void SetUsedFlag();

    void Merge(Contig &contig);
    bool Append(const ContigPath &path)
    {
        if (nodes.back() != path.nodes.front())
            return false;
        for (unsigned i = 1; i < path.nodes.size(); ++i)
            Append(path.nodes[i]);
        return true;
    }

    bool GetSuffix(const ContigNodeAdapter &begin, ContigPath &path)
    {
        path.Clear();
        unsigned index = 0;
        while (index < nodes.size() && nodes[index] != begin)
            ++index;

        if (index == nodes.size())
            return false;

        for (unsigned i = index; i < nodes.size(); ++i)
            path.Append(nodes[i]);

        return true;
    }

    bool Contain(ContigPath &path)
    {
        unsigned index = 0;
        while (index < nodes.size() && nodes[index] != path.nodes[0])
            ++index;

        if (index == nodes.size())
            return false;

        if (nodes.size() - index < path.nodes.size())
            return false;

        for (unsigned i = 0; i < path.nodes.size(); ++i)
        {
            if (nodes[i+index] != path.nodes[i])
                return false;
        }

        return true;
    }

protected:
    static const uint64 UnitOne = 100000000ULL;

    std::vector<ContigNodeAdapter> nodes;
    uint64 weight;
};

class GappedContigPath: public ContigPath
{
public:
    void Append(const ContigNodeAdapter &node, int d)
    { ContigPath::Append(node); if (nodes.size() > 1) distances.push_back(d); }

    void Pop()
    { ContigPath::Pop(); distances.pop_back(); }

    void ReverseComplement()
    { ContigPath::ReverseComplement(); std::reverse(distances.begin(), distances.end()); }

    void Merge(Contig &contig);

protected:
    std::vector<int> distances;
};

class ContigBranchGroup
{
public:
    ContigBranchGroup(ContigGraph *graph, ContigNodeAdapter &begin, int max_branches = 2, int max_length = kmerLength + 2)
    {
        this->graph = graph;
        this->begin = begin;
        this->max_branches = max_branches;
        this->max_length = max_length;
    }

    bool Search();
    void Merge();

private:
    ContigGraph *graph;
    ContigNodeAdapter begin;
    std::vector<ContigPath> branches;
    int max_branches;
    int max_length;
};

#endif
