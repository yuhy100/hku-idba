#include "globals.h"
#include "ScaffoldGraph.h"
#include "ContigBranchGroup.h"

#include <vector>
#include <cmath>

using namespace std;

static bool Compare(const ContigNode &x, const ContigNode &y)
{
    return x.GetSize() < y.GetSize();
}

void ScaffoldGraph::Initialize(vector<Contig> &contigs)
{
    fill_n(hits[0], MaxDistance, 0);
    fill_n(hits[1], MaxDistance, 0);

    sort(contigs.begin(), contigs.end(), Compare);
    reverse(contigs.begin(), contigs.end());

    ContigGraph::Initialize(contigs);
    BuildAllEdges();
    in_connections.resize(nodes.size());
    out_connections.resize(nodes.size());
    in_paths.resize(nodes.size());
    out_paths.resize(nodes.size());

    prev = prev_table + kmerLength;
}

void ScaffoldGraph::AddPair(vector<Alignment> &alignments1, vector<Alignment> &alignments2, int type)
{
    int index = 0;
    for (unsigned i = 0; i < alignments1.size(); ++i)
    {
        if (alignments1[i].contigLength < 2 * kmerLength)
            continue;

        alignments1[index++] = alignments1[i];
    }
    alignments1.resize(index);

    index = 0;
    for (unsigned i = 0; i < alignments2.size(); ++i)
    {
        if (alignments2[i].contigLength < 2 * kmerLength)
            continue;

        alignments2[index++] = alignments2[i];
    }
    alignments2.resize(index);


    //cout << alignments1.size() << " " << alignments2.size() << endl;
    if (alignments1.size() >= 2 && alignments2.size() >= 2)
        return;

    if (alignments1.size() > 1 || alignments2.size() > 1)
        return;

    for (unsigned i = 0; i < alignments1.size(); ++i)
    {
        for (unsigned j = 0; j < alignments2.size(); ++j)
        {
            AddPair(alignments1[i], alignments2[j], type);
        }
    }
}

void ScaffoldGraph::AddPair(Alignment a1, Alignment a2, int type)
{
    a2.ReverseComplement();

//    cout << a1.contigId << " " << a1.isReverse << " "
//        << a2.contigId << " " << as2.isReverse << endl;

    if (a1.contigId == a2.contigId && a1.isReverse == a2.isReverse)
    {
        int from = a1.contigOffset - a1.readOffset;
        int to = a2.contigOffset + (a2.readLength - a2.readOffset);
        if (!(to - from < 0 || to - from >= MaxDistance))
        {
            ++hits[type][to - from];
        }
    }
    else
    {
        int d = - (a1.contigLength - (a1.contigOffset - a1.readOffset))
                - (a2.contigOffset + (a2.readLength - a2.readOffset));

        AddConnection(ContigNodeAdapter(&nodes[a1.contigId], a1.isReverse),
                ContigNodeAdapter(&nodes[a2.contigId], a2.isReverse), d, type);
    }
}

void ScaffoldGraph::ComputeDistance(int type)
{
    int64 count = 0;
    for (int i = 0; i < MaxDistance; ++i)
        count += hits[type][i];

    int discard = count/200;
    int from = 0;
    int sum = 0;
    while (sum + hits[type][from] < discard)
    {
        sum += hits[type][from];
        ++from;
    }

    int to = MaxDistance - 1;
    sum = 0;
    while (sum + hits[type][to] < discard)
    {
        sum += hits[type][to];
        --to;
    }

    double sum_distance = 0;
    int real_num = 0;
    for (int i = from; i <= to; ++i)
    {
        sum_distance += i * hits[type][i];
        real_num += hits[type][i];
    }

    estimate_distance[type] = int(round(sum_distance/real_num));

    int region = real_num * 90 / 100;
    sum = hits[type][estimate_distance[type]];
    int offset = 1;
    while (offset < estimate_distance[type] && sum  < region)
    {
        sum += hits[type][estimate_distance[type] - offset] + hits[type][estimate_distance[type] + offset];
        ++offset;
    }

    if (offset < estimate_distance[type]/20)
        offset = estimate_distance[type]/20;

    delta[type] = offset;

    cout << estimate_distance[type] << " " << delta[type] << endl;
}

void ScaffoldGraph::FindUniquePaths()
{
    cout << "process connection" << endl;
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        ProcessConnections(in_connections[i]);
        ProcessConnections(out_connections[i]);
    }
    cout << "process connection" << endl;

    int conn_count = 0;
    int unique_count = 0;
    int unique_multi_degree = 0;
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        ContigNodeAdapter adp(&nodes[i]);
        for (int strand = 0; strand < 2; ++strand)
        {
            if (GetConnections(adp).size() > 0 && IsConsistance(adp, GetConnections(adp).back().to))
                    //&& GetConnections(adp).back().to.GetSize() > estimate_distance)
            {
                ++conn_count;
                ContigPath path;

                if (FindPath(adp, path))
                {
                    GetPaths(adp).push_back(path);
                    ++unique_count;

                    if (GetConnections(adp).size() > 1)
                        ++unique_multi_degree;
                }

            }

            adp.ReverseComplement();
        }
    }

    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        ProcessPaths(in_paths[i]);
        ProcessPaths(out_paths[i]);
    }

    cout << unique_count << " " << conn_count << " " << unique_multi_degree << endl;
}

int64 ScaffoldGraph::Scaffold(vector<Contig> &contigs)
{
    contigs.resize(0);

    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i].IsDead() || nodes[i].IsUsed())
            continue;

        ContigNodeAdapter contig_adp(&nodes[i]);
        contig_adp.SetUsedFlag();
        ContigPath contig_path;
        contig_path.Append(contig_adp);
        for (int strand = 0; strand < 2; ++strand)
        {
            while (true)
            {
                ContigNodeAdapter curr = contig_path.GetEndNodeAdapter();
                ContigPath path;

                if (GetPaths(curr).size() != 1)
                    break;

                path = GetPaths(curr)[0];

                while (true)
                {
                    ContigNodeAdapter end = path.GetEndNodeAdapter();
                    ContigNodeAdapter tmp_begin = end;
                    tmp_begin.ReverseComplement();

                    if (path.Size() > 1 &&
                            (end.IsUsed() || GetPaths(tmp_begin).size() != 1))
                        path.Pop();
                    else
                        break;
                }

                ContigPath tmp_path;
                ContigNodeAdapter tmp_end = curr;
                ContigNodeAdapter tmp_begin = path.GetEndNodeAdapter();
                tmp_begin.ReverseComplement();
                tmp_end.ReverseComplement();

                if (GetPaths(tmp_begin).size() != 1)
                    break;

                if (path.GetEndNodeAdapter().IsUsed())
                    break;

                path.SetUsedFlag();
                contig_path.Append(path);
            }

            contig_path.ReverseComplement();
        }

        Contig contig;
        contig_path.Merge(contig);
        contigs.push_back(contig);
    }

    return contigs.size();
}

int64 ScaffoldGraph::ScaffoldWithGap(vector<Contig> &contigs)
{
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        ProcessLongConnections(in_connections[i]);
        ProcessLongConnections(out_connections[i]);
    }

    contigs.resize(0);

    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i].IsDead() || nodes[i].IsUsed())
            continue;

        ContigNodeAdapter contig_adp(&nodes[i]);

//        if (!IsLongSeed(contig_adp))
//            continue;

        contig_adp.SetUsedFlag();
        GappedContigPath contig_path;
        contig_path.Append(contig_adp, 0);
        for (int strand = 0; strand < 2; ++strand)
        {
            while (true)
            {
                ContigNodeAdapter curr = contig_path.GetEndNodeAdapter();

                if (GetConnections(curr).size() != 1)
                    break;

                ContigNodeAdapter next = GetConnections(curr)[0].to;
                ContigNodeAdapter temp_begin = next;
                temp_begin.ReverseComplement();

                if (GetConnections(temp_begin).size() != 1)
                    break;

                if (next.IsUsed())
                    break;

                next.SetUsedFlag();
                contig_path.Append(next, GetConnections(curr)[0].distance);
            }

            contig_path.ReverseComplement();
        }

        Contig contig;
        contig_path.Merge(contig);
        contigs.push_back(contig);
    }

    return contigs.size();
}

void ScaffoldGraph::ProcessConnections(vector<ScaffoldConnection> &connections)
{
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        connections[i].Translate(estimate_distance[0], 0);
        connections[i].Translate(estimate_distance[1], 1);
        connections[i].ComputeDistance();
    }

    int index = 0;
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        if ((int)connections[i].values.size() >= min_pairs 
                && connections[i].distance >= -2*kmerLength 
                //&& connections[i].distance <= estimate_distance[connections[i].type] * 2
                && connections[i].IsConsistant(delta[0]))
            connections[index++] = connections[i];
    }
    connections.resize(index);

    sort(connections.begin(), connections.end());
}

void ScaffoldGraph::ProcessPaths(std::vector<ContigPath> &paths)
{
    int index = 0;
    for (unsigned i = 0; i < paths.size(); ++i)
    {
        bool flag = true;
        for (unsigned j = i+1; j < paths.size(); ++j)
        {
            if (paths[j].Contain(paths[i]))
            {
                flag = false;
                break;
            }
        }

        if (flag)
            paths[index++] = paths[i];
    }
    paths.resize(index);
}

void ScaffoldGraph::ProcessLongConnections(vector<ScaffoldConnection> &connections)
{
    int index = 0;
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        if (IsLongSeed(connections[i].to))
            connections[index++] = connections[i];
    }
    connections.resize(index);
}

bool ScaffoldGraph::FindPath(ContigNodeAdapter &node, ContigPath &path)
{
    if (GetConnections(node).size() == 0)
        return false;

    if (!IsConsistance(node, GetConnections(node).back().to))
        return false;

    //return FindPath(node, GetConnections(node).back(), path, GetConnections(node).back().distance);

    path.Append(node);

    vector<ScaffoldConnection> &connections = GetConnections(node);
    int length = 0;
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        ContigPath middle_path;
        if (!FindPath(path.GetEndNodeAdapter(), connections[i].to, middle_path, connections[i].distance - length))
            return false;

        length = connections[i].distance + connections[i].to.GetSize();

        path.Append(middle_path);
    }

    return true;
}

bool ScaffoldGraph::FindPath(ContigNodeAdapter &node, ContigNodeAdapter &target, ContigPath &path, int length)
{
    for (int i = -kmerLength; i <= length + delta[0]; ++i)
        prev[i].clear();

    prev[-kmerLength][node].push_back(node);

    static int index = 0;
    cout << index++ << " " << node.Data() << " " << target.Data() << " " << node.GetSize() << " " << target.GetSize() << " " << length << endl;

    used_time = 0;
    for (int i = -kmerLength; i <= length + delta[0]; ++i)
    {
        if (i >= length)
        {
            int offset = i - length;
            if (length - offset >= -kmerLength && prev[length - offset].find(target) != prev[length - offset].end())
                break;

            if (prev[i].find(target) != prev[i].end())
                break;
        }

//        if (length == 591)
//            cout << i << " " << prev[i].size() << endl;
        used_time += prev[i].size();
        if (used_time >= TimeLimit)
        {
            cout << "time out " << length << endl;

//            for (int d = -kmerLength; d <= length + delta[0]; ++d)
//            {
//                cout << "distance " << d << endl;
//
//                map<ContigNodeAdapter, vector<ContigNodeAdapter> >::iterator iter = prev[d].begin();
//                while (iter != prev[d].end())
//                {
//                    ContigNodeAdapter current = iter->first;
//                    int dist;
//                    if (current == node)
//                        dist = - kmerLength + 1;
//                    else
//                        dist = d - kmerLength + 1 + current.GetSize();
//
//                    Contig contig;
//                    current.GetContig(contig);
//
//                    contig.Decode();
//                    cout << current.Data() << " " << dist << " " << current.GetSize() << " " << contig << endl;
//
//                    ++iter;
//                }
//
//                cout << endl;
//            }

            break;
        }

        map<ContigNodeAdapter, vector<ContigNodeAdapter> >::iterator iter = prev[i].begin();
        while (iter != prev[i].end())
        {
            ContigNodeAdapter current = iter->first;
            vector<ContigNodeAdapter> neighbors;
            GetNeighbors(current, neighbors);

            int d;
            if (current == node)
                d = - kmerLength + 1;
            else
                d = i - kmerLength + 1 + current.GetSize();

            if (d <= length + delta[0])
            {
                for (unsigned j = 0; j < neighbors.size(); ++j)
                    prev[d][neighbors[j]].push_back(current);
            }

            ++iter;
        }
    }

    cout << "dp end " << length << endl;


    int result_length = MaxDistance;

    for (int i = 0; i <= delta[0]; ++i)
    {
        if (length + i >= -kmerLength && prev[length + i].find(target) != prev[length + i].end())
        {
            result_length = length + i;
            break;
        }

        if (length - i >= -kmerLength && prev[length - i].find(target) != prev[length - i].end())
        {
            result_length = length - i;
            break;
        }
    }

    cout << "result length " <<  result_length << endl;

    if (result_length == MaxDistance)
        return false;

//    cout << node.GetSize() << endl;
//    cout << "true" << endl;
    result_paths.resize(0);
    if (!ConstructPath(target, result_length))
        return false;
    path = result_paths[0];

    static int long_count = 0;
    if (length > 200)
    {
        cout << "long " << long_count++ << " " << length << endl;
    }

    return true;
}

bool ScaffoldGraph::ConstructPath(ContigNodeAdapter &target, int length)
{
    partial_path.push_back(target);

    if (length == -target.GetSize() && length <= -kmerLength)
    {
        //cout << "good" << endl;
        reverse(partial_path.begin(), partial_path.end());
        current_path.Clear();
        for (unsigned i = 0; i < partial_path.size(); ++i)
        {
            current_path.Append(partial_path[i]);
        }
        result_paths.push_back(current_path);
        reverse(partial_path.begin(), partial_path.end());

        partial_path.pop_back();

        if (result_paths.size() > 10)
        {
            cout << "too many path" << endl;
            return false;
        }

        return true;
    }

    vector<ContigNodeAdapter> &v = prev[length][target];
    for (unsigned i = 0; i < v.size(); ++i)
    {
        int d = length - (-kmerLength + 1 + v[i].GetSize());
        if (!ConstructPath(v[i], d))
        {
            partial_path.pop_back();
            return false;
        }
    }

    partial_path.pop_back();

    return true;
}

bool ScaffoldGraph::FindPath(ContigNodeAdapter &node, ScaffoldConnection &connection, ContigPath &path, int length)
{
    result_paths.resize(0);
    used_time = 0;
    expected_length = length;

    DepthFirstSearch(node, connection.to);

    if (used_time >= TimeLimit)
        return false;

    if (result_paths.size() == 1)
    {
        int best = 0;
        for (unsigned i = 0; i < result_paths.size(); ++i)
        {
            if (result_paths[i].Weight() > result_paths[best].Weight())
                best = i;
        }
        path = result_paths[best];
        //cout << expected_length << " " << path.InternalDistance() << endl;
        return true;
    }
    else if (result_paths.size() > 1)
    {
        int index = 0;
        for (unsigned i = 0; i < result_paths.size(); ++i)
        {
            Contig contig;
            result_paths[i].Merge(contig);
            contig.Decode();

            bool flag = true;
            for (unsigned j = i+1; j < result_paths.size(); ++j)
            {
                Contig contig2;
                result_paths[j].Merge(contig2);
                contig2.Decode();

                int size = min(contig.Size(), contig2.Size());
                int count = 0;
                for (int i = 0; i < size; ++i)
                {
                    if (contig[i] != contig2[i])
                        ++count;
                }

                if (abs(contig.Size() - contig2.Size()) < 0.01 * size
                        && count < 0.01 * size)
                    flag = false;
            }

            if (flag)
                result_paths[index++] = result_paths[i];
        }

        result_paths.resize(index);

        int best = 0;
        int minimum = 1000000;
        for (unsigned i = 0; i < result_paths.size(); ++i)
        {
            int d = abs(result_paths[i].InternalDistance() - expected_length);
            if (d < minimum)
            {
                best = i;
                minimum = d;
            }
        }

        path = result_paths[best];
        return true;
    }
    else
        return false;
}

void ScaffoldGraph::DepthFirstSearch(ContigNodeAdapter &node, ContigNodeAdapter &target)
{
    if (node == target)
    {
        current_path.Append(node);
        if (abs(current_path.InternalDistance() - expected_length) <= delta[0])
        {
            //cout << current_path.InternalDistance() << " " << expected_length << endl;
            result_paths.push_back(current_path);
        }
        current_path.Pop();

        return;
    }

    ++used_time;
    if (used_time >= TimeLimit)
        return;

    if (current_path.InternalDistance() > expected_length + delta[0])
        return;

    vector<ContigNodeAdapter> neighbors;
    GetNeighbors(node, neighbors);

    current_path.Append(node);
    node.SetDeadFlag();

    for (unsigned i = 0; i < neighbors.size(); ++i)
    {
        //if (!neighbors[i].IsDead())
        {
            DepthFirstSearch(neighbors[i], target);
        }
    }

    node.ResetDeadFlag();
    current_path.Pop();
}

bool ScaffoldGraph::IsConsistance(const ContigNodeAdapter &node, const ContigNodeAdapter &end)
{
    return true;
    vector<ScaffoldConnection> &connections = GetConnections(node);
    if (connections.back().to != end)
        return false;

    for (unsigned i = 0; i+1 < connections.size(); ++i)
    {
        ContigNodeAdapter middle = connections[i].to;
        ContigNodeAdapter next = connections[i+1].to;

        int d = GetDistance(node, middle) + node.GetSize() + GetDistance(middle, next);

        if (abs(d - GetDistance(node, next) > delta[0]))
        {
            //cout << d << " " << GetDistance(node, next) << endl;
            return false;
        }
    }

    return true;
}


