#include "globals.h"
#include "kmer.h"
#include "sequence.h"
#include "hashTable.h"
#include "utils.h"
#include "contigGraph.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

// bool compareContig(int x, int y)
// {
//     return contigs[x].length > contigs[y].length;
// }

// static int GetOverlap(Sequence *s1, Sequence *s2)
// {
//     for (int i = kmerLength - 1; i > 10; --i)
//     {
//         Kmer a = s1->GetKmer(0, i);
//         Kmer b = s2->GetKmer(s2->length - i, i);
//         if (a == b)
//             return i;
//     }
//
//     return 0;
// }

struct CompareContig
{
    ContigGraph *contigGraph;
    CompareContig(ContigGraph *contigGraph) { this->contigGraph = contigGraph; }
    bool operator() (int x, int y)
    {
        return contigGraph->contigs[x].length > contigGraph->contigs[y].length;
    }
};

static bool comparePath(const Path &p1, const Path &p2)
{
    return p1.points.size() < p2.points.size();
}

void ContigGraph::ReadFromFile(FILE *fp)
{
    numContigs = CountSequences(fp);
    contigs = new Sequence[numContigs*2];
    connections = new vector<Connection>[numContigs*2];
    paths = new vector<Path>[numContigs*2];
    possible = new vector<Path>[numContigs*2];
    target = new int[numContigs*2];
    fill_n(target, numContigs*2, 0);

    Sequence seq;
    for (unsigned i = 0; i < numContigs; ++i)
    {
        ReadFasta(fp, &seq);
        seq.Encode();
        contigs[i<<1] = seq;
        seq.ReverseComplement();
        contigs[(i<<1)|1] = seq;
    }
}

void ContigGraph::AddConnection(int from, int to, int d)
{
    bool exist = false;
    for (unsigned i = 0; i < connections[from].size(); ++i)
    {
        if (connections[from][i].to == to)
        {
            connections[from][i].values.push_back(d);
            exist = true;
            break;
        }
    }

    if (!exist)
    {
        Connection connection;
        connection.from = from;
        connection.to = to;
        connection.values.push_back(d);
        connection.distance = 0;
        connections[from].push_back(connection);
    }
}

void ContigGraph::FilterConnections(int minPairs)
{
//     int count = 0;
//     int nCount = 0;
//     int unique = 0;
    for (unsigned x = 0; x < numContigs*2; ++x)
    {
//         int branch = 0;
        int l = 0;
        for (unsigned j = 0; j < connections[x].size(); ++j)
        {
            if (connections[x][j].values.size() >= (unsigned)minPairs)// && contigs[x].length > 100 && contigs[edges[x][j].to].length > 100)
            {
                int sum = 0;
                sort(connections[x][j].values.begin(), connections[x][j].values.end());
                for (unsigned k = 0; k < connections[x][j].values.size(); ++k)
                {
                    sum += connections[x][j].values[k];
                }
                double d = double(sum)/(connections[x][j].values.size());

                bool flag = true;
                for (unsigned k = 0; k < connections[x][j].values.size(); ++k)
                {
                    if (fabs(d - connections[x][j].values[k]) > range)
                    {
                        flag = false;
                        break;
                    }
                }

                if (!flag)
                    continue;

//                 ++count;
//                 if (d < 0)
//                     ++nCount;
//                 ++branch;

                connections[x][j].distance = int (round(d));
                connections[x][l++] = connections[x][j];
            }
        }
        connections[x].resize(l);

//         if (branch == 1)
//             ++unique;
    }

//     cout << count << " " << nCount << " " << unique << endl;
}

void ContigGraph::BuildInitialPaths()
{
    HashTable *hashTable = new HashTable(numContigs << 2);
    vector<vector<int> > positions;

    kmerLength--;
    for (unsigned i = 0; i < numContigs*2; ++i)
    {
        Kmer kmer = contigs[i].GetKmer(0, kmerLength);

        HashNode *p = hashTable->InsertKmer(kmer);
        if (p->count == 1)
        {
            p->Data() = hashTable->numNodes - 1;
            positions.resize(positions.size() + 1);
        }
        positions[p->Data()].push_back(i);
    }

    for (unsigned i = 0; i < numContigs*2; ++i)
    {
        Kmer kmer = contigs[i].GetKmer(contigs[i].length - kmerLength, kmerLength);
        HashNode *p = hashTable->GetByKmer(kmer);

        if (p == NULL)
            continue;

        vector<int> &v = positions[p->Data()];
        for (unsigned k = 0; k < v.size(); ++k)
        {
            unsigned j = v[k];
            if ((i >> 1) == (j >> 1))
                continue;

//             Kmer b = contigs[j].GetKmer(0, kmerLength);
//             if (a == b)
            {
                Path path;
                path.points.push_back(j);
                path.from = i;
                path.to = j;
                path.length = -kmerLength;
                paths[i].push_back(path);
            }
        }
    }
    kmerLength++;
}

void ContigGraph::DepthFirstSearch(int from, Path &path)
{
    if (--timeLimit < 0)
        return ;

    if (path.length > pathLimit)
        return;

    if (target[from] == 1)
    {
        target[from] = 2;
        ++found;
    }

    if (found == totalTarget)
    {
        result.push_back(path);

        if (target[from] == 2)
        {
            target[from] = 1;
            --found;
        }

        return;
    }

    for (unsigned i = 0; i < paths[from].size(); ++i)
    {
        Path &p = paths[from][i];
        for (unsigned j = 0; j < p.points.size(); ++j)
            path.points.push_back(p.points[j]);
        path.length += contigs[from].length + p.length;
        path.to = p.to;

        DepthFirstSearch(p.to, path);

        for (unsigned j = 0; j < p.points.size(); ++j)
            path.points.pop_back();
        path.length -= contigs[from].length + p.length;
    }

    if (target[from] == 2)
    {
        target[from] = 1;
        --found;
    }
}

void ContigGraph::FindPaths(int x, std::vector<Path> &outPaths)
{
    int maximum = -1000000;
//     int multi = 0;
    for (unsigned i = 0; i < connections[x].size(); ++i)
    {
        target[connections[x][i].to] = 1;
        if (connections[x][i].distance > maximum)
            maximum = connections[x][i].distance;
    }

    totalTarget = connections[x].size();
    found = 0;
    result.resize(0);
    pathLimit = maximum + range;
    timeLimit = TimeLimit;
    Path path;
    path.from = x;
    path.to = x;
    path.length = -contigs[x].length;

    DepthFirstSearch(x, path);

//     if (result.size() > 0)
//     {
//         ++multi;
//     }

    if (result.size() > 0)
    {
        int l = 0;
        for (unsigned i = 0; i < result.size(); ++i)
        {
            int length = -contigs[x].length;
            int k = x;
            for (unsigned j = 0; j < result[i].points.size(); ++j)
            {
                int y = result[i].points[j];
                length += contigs[k].length - (kmerLength - 1);
                if (target[y])
                    target[y] = length;
                k = y;
            }

            bool flag = true;
            for (unsigned j = 0; j < connections[x].size(); ++j)
            {
                int y = connections[x][j].to;
                if (!(target[y] >= connections[x][j].distance - range
                      && target[y] <= connections[x][j].distance + range))
                {
                    flag = false;
                    break;
                }
            }

            if (flag)
                result[l++] = result[i];
        }
        result.resize(l);
    }

    for (unsigned i = 0; i < connections[x].size(); ++i)
    {
        target[connections[x][i].to] = 0;
    }

    int l = 0;
    for (unsigned i = 0; i < result.size(); ++i)
    {
        bool flag = true;
        for (int j = 0; j < l; ++j)
        {
            Sequence seq1;
            Sequence seq2;
            for (unsigned k = 0; k+1 < result[i].points.size(); ++k)
            {
                int y = result[i].points[k];
                for (int a = kmerLength-1; a < contigs[y].length; ++a)
                    seq1 += contigs[y].bases[a];
            }
            for (unsigned k = 0; k+1 < result[j].points.size(); ++k)
            {
                int y = result[j].points[k];
                for (int a = kmerLength-1; a < contigs[y].length; ++a)
                    seq2 += contigs[y].bases[a];
            }

            int diff = 0;
            for (int k = 0; k < min(seq1.length, seq2.length); ++k)
            {
                if (seq1[k] != seq2[k])
                    ++diff;
            }

            int maximum = max(seq1.length, seq2.length);
            if (abs(seq1.length - seq2.length) < maximum * 0.01 && diff < 0.01)
            {
                flag = false;
                break;
            }
        }

        if (flag)
            result[l++] = result[i];
    }
    result.resize(l);

    if (timeLimit < 0)
        result.resize(0);

    outPaths.swap(result);
//     outPaths = result;
}

void ContigGraph::FindPossibleConnections()
{
//     int multi = 0;
    for (unsigned x = 0; x < numContigs*2; ++x)
    {
        if (connections[x].size() == 0 && contigs[x].length < 2*kmerLength)
            continue;

        vector<Path> outPaths;
        FindPaths(x, outPaths);

        if (outPaths.size() == 1)
        {
            possible[x].push_back(outPaths[0]);

            vector<int> v = outPaths[0].points;
            reverse(v.begin(), v.end());
            v.push_back(x);
            for (unsigned i = 0; i < v.size(); ++i)
                v[i] ^= 1;

            for (unsigned i = 0; i+1 < v.size(); ++i)
            {
                Path p;
                p.from = v[i];
                p.to = v.back();
                p.length = 0;
                for (unsigned j = i+1; j < v.size(); ++j)
                {
                    p.points.push_back(v[j]);
                    p.length += contigs[v[j]].length - (kmerLength - 1);
                }
                p.length -= contigs[v.back()].length;

                if (contigs[p.from].length > kmerLength*2)
                    possible[p.from].push_back(p);
            }
        }
        else
        {
//             if (connections[x].size() == 1 && connections[x][0].distance < 0)
//             {
//                 int y = connections[x][0].to;
//                 int d = connections[x][0].distance;
//
//                 int overlap = GetOverlap(&contigs[x], &contigs[y]);
//                 if (overlap > 20 && overlap >= -d - range && overlap <= -d + range)
//                 {
//                     Path path;
//                     path.from = x;
//                     path.to = y;
//                     path.length = -overlap;
//                     path.points.push_back(y);
//                     possible[x].push_back(path);
//                 }
//             }
        }
    }

    int pairs = 0;
    vector<int> v;
    for (unsigned i = 0; i < numContigs*2; i += 2)
        v.push_back(i);
//     sort(v.begin(), v.end(), CompareContig(this));

    for (int i = v.size()-1; i >= 0; --i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            int x = v[i] + strand;

            if (possible[x].size() == 0)
                continue;

            sort(possible[x].begin(), possible[x].end(), comparePath);

            bool flag = true;
            for (unsigned i = 0; i+1 < possible[x].size(); ++i)
            {
                for (unsigned j = 0; j < possible[x][i].points.size(); ++j)
                {
                    if (possible[x][i].points[j] != possible[x][i+1].points[j])
                    {
                        flag = false;
                        break;
                    }
                }

                if (!flag)
                    break;
            }

            if (flag)
            {
                possible[x][0] = possible[x].back();
                possible[x].resize(1);
            }

            if (possible[x].size() == 1)
            {
                ++pairs;
                paths[x].resize(1);
                paths[x][0] = possible[x][0];
            }
            else if (possible[x].size() > 1)
            {
                paths[x] = possible[x];
            }
        }
    }

//     cerr << "pairs " << pairs << endl;
}

void ContigGraph::MergeContigs(std::vector<Sequence> &mergedContigs)
{

    char *isUsed = new char[numContigs*2];
    fill_n(isUsed, numContigs*2, 0);

//     for (unsigned k = 0; k < v.size(); ++k)
    for (unsigned i = 0; i < numContigs*2; ++i)
    {
//         int i = v[k];
        if (isUsed[i>>1])
            continue;

//         if (contigs[i].length < 2*kmerLength)// || contigs[i].length < 100)
//             continue;

        Sequence contig = contigs[i];
        isUsed[i>>1] = 1;
        for (int strand = 0; strand < 2; ++strand)
        {
            int x = i + strand;
            while (paths[x].size() == 1)
            {
                Path &path = paths[x][0];

                if (path.points.size() > 1)
                {
                    while (path.points.size() > 1 && (paths[path.to^1].size() != 1 || possible[path.to^1].size() > 1 || isUsed[path.to>>1]))
                    {
                        path.points.pop_back();
                        path.length -= (contigs[path.points.back()].length - (kmerLength -1));
                        path.to = path.points.back();
                    }

                    if (isUsed[path.to>>1])
                        break;

                    if (paths[path.to^1].size() != 1 || possible[path.to^1].size() > 1)
                        break;

                    for (unsigned j = 0; j+1 < path.points.size(); ++j)
                    {
                        int y = path.points[j];
                        isUsed[y>>1] = 1;
                        for (int k = kmerLength-1; k < contigs[y].length; ++k)
                            contig += contigs[y].bases[k];
                    }

                    if (paths[path.to^1].size() != 1 || possible[path.to^1].size() > 1)
                        break;

                    if (isUsed[path.to>>1])
                        break;

                    int y = path.to;
                    isUsed[y>>1] = 1;
                    for (int k = kmerLength-1; k < contigs[y].length; ++k)
                        contig += contigs[y].bases[k];
                }
                else
                {
                    if (isUsed[path.to>>1])
                        break;

                    if (paths[path.to^1].size() != 1 || possible[path.to^1].size() > 1)
                        break;

                    int y = path.to;
                    isUsed[y>>1] = 1;
                    for (int k = -path.length; k < contigs[y].length; ++k)
                        contig += contigs[y].bases[k];
                }

                x = path.to;
                isUsed[x>>1] = 1;
            }

            contig.ReverseComplement();
        }

        contig.Decode();
        mergedContigs.push_back(contig);
    }
}
