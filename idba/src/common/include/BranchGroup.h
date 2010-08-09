#ifndef __BRANCH_GROUP_H_

#define __BRANCH_GROUP_H_

#include "globals.h"
#include "Kmer.h"
#include "HashNode.h"

#include <vector>

class HashGraph;

class Path
{
public:
    Path() { weight = 0; }

    void Append(KmerNodeAdapter &next)
    {
        path.push_back(next);
        weight += UnitOne / next.Count();
    }

    KmerNodeAdapter GetEndNodeAdapter() { return path.back(); }
    int Size() { return path.size(); }
    int64 Weight() { return weight; }

    bool IsSimplePath();
    void Inactivate();
    void Activate();

private:
    static const uint64 UnitOne = 100000000ULL;

    std::vector<KmerNodeAdapter> path;
    uint64 weight;
};

class BranchGroup
{
public:
    BranchGroup(HashGraph *graph, KmerNodeAdapter &begin, int max_branches = 2, int max_length = kmerLength + 2)
    {
        this->graph = graph;
        this->begin = begin;
        this->max_branches = max_branches;
        this->max_length = max_length;
    }

    bool Search();
    void Merge();


private:
    HashGraph *graph;
    KmerNodeAdapter begin;
    std::vector<Path> branches;
    int max_branches;
    int max_length;
};

#endif
