#ifndef __MAP_GRAPH_H_

#define __MAP_GRAPH_H_

#include "globals.h"
#include "Contig.h"

#include <cstdio>
#include <vector>

struct Sequence;

struct Connection
{
    int from;
    int to;
    int distance;
    std::vector<int> values;
};

struct Path
{
    int from;
    int to;
    int length;
    std::vector<int> points;
};

struct MapGraph
{
    Sequence *contigs;
    std::vector<Connection> *connections;
    std::vector<Path> *paths;
    std::vector<Path> *possible;
    std::vector<Path> result;
    uint64 numContigs;
    int *target;
    int totalTarget;
    int timeLimit;
    int pathLimit;
    int found;
    int range;
    int distance;

    MapGraph()
    {
        contigs = NULL;
        connections = NULL;
        paths = NULL;
        possible = NULL;
        target = NULL;
        range = 0;
    }

    ~MapGraph()
    {
        delete [] contigs;
        delete [] connections;
        delete [] paths;
        delete [] possible;
        delete [] target;
    }

    bool ValidatePath(int x, Connection &connection);
    void ReadFromFile(std::FILE *fp);
    void Initialize(std::vector<Contig> &contigs);
    void SetRange(int range) { this->range = range; }
    void SetDistance(int distance) { this->distance = distance; }

    void AddConnection(int from, int to, int d);
    void FilterConnections(int minPairs);
    void BuildInitialPaths();
    void DepthFirstSearch(int from, Path &path);
    void FindPaths(int from, std::vector<Path> &outPaths);
    void FindPossibleConnections();
    void MergeContigs(std::vector<Sequence> &mergedContigs);

    static const uint64 TimeLimit = 10000;
};

#endif
