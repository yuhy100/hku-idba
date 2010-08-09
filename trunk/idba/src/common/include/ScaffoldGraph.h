#ifndef __SCAFFOLD_GRAPH_H_

#define __SCAFFOLD_GRAPH_H_

#include "globals.h"
#include "ContigNode.h"
#include "ContigGraph.h"
#include "ContigBranchGroup.h"
#include "HashAlign.h"

#include <algorithm>
#include <vector>
#include <iostream>
#include <map>

class ScaffoldConnection
{
public:
    ContigNodeAdapter to;
    int distance;
    std::vector<int> values;
    std::vector<int> types;

    ScaffoldConnection() { }

    ScaffoldConnection(const ContigNodeAdapter &node, int value, int type = 0)
    { to = node; values.push_back(value); types.push_back(type); }

    bool operator <(const ScaffoldConnection &connection) const
    { return distance < connection.distance; }

    bool AddValue(const ContigNodeAdapter &node, int value, int type = 0)
    {
        if (node != to)
            return false;
        values.push_back(value);
        types.push_back(type);
        return true;
    }

    void ComputeDistance()
    {
        std::sort(values.begin(), values.end());
        distance = values[values.size()/2];
    }

    void Translate(int offset, int type = 0)
    {
        for (unsigned i = 0; i < values.size(); ++i)
        {
            if (types[i] == type)
                values[i] += offset;
        }
    }

    bool IsConsistant(int delta)
    {
       int ignore = values.size()/10;
       for (unsigned i = ignore; i + ignore < values.size(); ++i)
       {
           if (abs(values[i] - distance) > delta)
               return false;
       }
       return true;
    }
};

class ScaffoldGraph: public ContigGraph
{
public:
    ScaffoldGraph(int min_pairs = 5) 
    { SetMinPairs(min_pairs); }
    ScaffoldGraph(std::vector<Contig> &contigs, int min_pairs = 5) 
    { Initialize(contigs); SetMinPairs(min_pairs); }
    ~ScaffoldGraph() {}

    void Initialize(std::vector<Contig> &contigs);

    void AddPair(std::vector<Alignment> &alignments1, std::vector<Alignment> &alignments2, int type = 0);
    void AddPair(Alignment a1, Alignment a2, int type = 0);
    void ComputeDistance(int type = 0);
    void FindUniquePaths();

    int64 Scaffold(std::vector<Contig> &contigs);
    int64 ScaffoldWithGap(std::vector<Contig> &contigs);

    void SetMinPairs(int min_pairs) 
    { this->min_pairs = min_pairs; }

private:
    ScaffoldGraph(const ScaffoldGraph &);
    const ScaffoldGraph &operator =(const ScaffoldGraph &);


    void AddConnection(ContigNodeAdapter node1, ContigNodeAdapter node2, int distance, int type)
    {
        if (!IsSeed(node1) || !IsSeed(node2))
            return;

        AddConnection(GetConnections(node1), node2, distance, type);
        node1.ReverseComplement();
        node2.ReverseComplement();
        AddConnection(GetConnections(node2), node1, distance, type);
        //GetConnections(node2).push_back(node1);
    }

    std::vector<ScaffoldConnection> &GetConnections(const ContigNodeAdapter &node)
    { return (!node.IsReverse() ? out_connections[node.Data()] : in_connections[node.Data()]); }

    const std::vector<ScaffoldConnection> &GetConnections(const ContigNodeAdapter &node) const
    { return (!node.IsReverse() ? out_connections[node.Data()] : in_connections[node.Data()]); }

    std::vector<ContigPath> &GetPaths(const ContigNodeAdapter &node)
    { return (!node.IsReverse() ? out_paths[node.Data()] : in_paths[node.Data()]); }

    const std::vector<ContigPath> &GetPaths(const ContigNodeAdapter &node) const
    { return (!node.IsReverse() ? out_paths[node.Data()] : in_paths[node.Data()]); }

    int GetDistance(const ContigNodeAdapter &node, const ContigNodeAdapter &next)
    {
        std::vector<ScaffoldConnection> &connections = GetConnections(node);
        for (unsigned i = 0; i < connections.size(); ++i)
        {
            if (connections[i].to == next)
                return connections[i].distance;
        }

        return MaxDistance;
    }

    void AddConnection(std::vector<ScaffoldConnection> &connections, 
            ContigNodeAdapter &node, int distance, int type)
    {
        for (unsigned i = 0; i < connections.size(); ++i)
        {
            if (connections[i].AddValue(node, distance, type))
                return;
        }

        connections.push_back(ScaffoldConnection(node, distance, type));
    }

    void ProcessConnections(std::vector<ScaffoldConnection> &connections);
    void ProcessPaths(std::vector<ContigPath> &paths);
    void ProcessLongConnections(std::vector<ScaffoldConnection> &connects);

    bool FindPath(ContigNodeAdapter &node, ContigPath &path);
    bool FindPath(ContigNodeAdapter &node, ContigNodeAdapter &target, ContigPath &path, int length);
    bool FindPath(ContigNodeAdapter &node, ScaffoldConnection &connection, ContigPath &path, int length);
    bool ConstructPath(ContigNodeAdapter &node, int length);
    void DepthFirstSearch(ContigNodeAdapter &node, ContigNodeAdapter &target);
    
    bool IsConsistance(const ContigNodeAdapter &node, const ContigNodeAdapter &end);
    bool IsSeed(const ContigNodeAdapter &node)
    { return node.GetSize() >= kmerLength * 2; }
    bool IsLongSeed(const ContigNodeAdapter &node)
    { return node.GetSize() >= estimate_distance[0] + delta[0]; }

    static const int MaxDistance = 10000;
    static const int TimeLimit = 10000;

    int expected_length;
    int used_time;
    std::map<ContigNodeAdapter, std::vector<ContigNodeAdapter> > *prev;
    std::map<ContigNodeAdapter, std::vector<ContigNodeAdapter> > prev_table[MaxDistance];
    ContigPath current_path;
    std::vector<ContigNodeAdapter> partial_path;
    std::vector<ContigPath> result_paths;
    std::vector<std::vector<ContigPath> > out_paths;
    std::vector<std::vector<ContigPath> > in_paths;
    std::vector<std::vector<ScaffoldConnection> > out_connections;
    std::vector<std::vector<ScaffoldConnection> > in_connections;

    int64 hits[2][MaxDistance];
    int estimate_distance[2];
    int delta[2];
    int min_pairs;
};

#endif
