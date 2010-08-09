#include "globals.h"
#include "Sequence.h"
#include "Utils.h"
#include "BitOperation.h"
#include "Reader.h"
#include "KmerVector.h"

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <map>


using namespace std;

const int MaxCluster = 2048;

int max_cluster_id = 0;
int flags[MaxCluster] = {0};
double distance[MaxCluster][MaxCluster];
vector<KmerVector> kmer_vectors[MaxCluster];
vector<int> childs[MaxCluster];
double min_distance[MaxCluster] = {0};

void Usage()
{
    fprintf(stderr, "usage: metaClusterHi group1.fa group2.fa ...\n");
}

double Distance(KmerVector *kv1, KmerVector *kv2)
{
    double sum = 0;
    for (int i = 0; i < kv1->size; ++i)
        sum += abs(kv1->ranks[i] - kv2->ranks[i]);
    return sum;
}

double DistanceMedian(vector<KmerVector> &kv1, vector<KmerVector> &kv2)
{
    vector<double> tmp(kv1.size() * kv2.size());
#pragma omp parallel for
    for (int i = 0; i < kv1.size(); ++i)
    {
        for (int j = 0; j < kv2.size(); ++j)
        {
            double d = Distance(&kv1[i], &kv2[j]);
            tmp[i*kv2.size() + j] = d;
        }
    }

    sort(tmp.begin(), tmp.end());
    return tmp[tmp.size()/2];
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        Usage();
        exit(1);
    }

    int num_cluster = argc - 1;
    for (int i = 0; i < num_cluster; ++i)
    {
        printf("%s\n", argv[i+1]);
        fflush(NULL);

        KmerVector kv;
        FastAReader reader(argv[i+1]);
        Sequence seq;
        string comment;

        reader.Read(seq, comment);

        while (reader.Read(seq, comment))
        {
            seq.Encode();
            kv.Compute(&seq);
            kv.ComputeRank();
            kmer_vectors[i].push_back(kv);
        }
    }

    int max_cluster_id = num_cluster;
    int num_round = num_cluster - 1;

    for (int round = 0; round < num_round; ++round)
    {
        int c1 = -1;
        int c2 = -1;
        double minimum = 1e100;

        for (int i = 0; i < max_cluster_id; ++i)
        {
            if (flags[i])
                continue;

            for (int j = i+1; j < max_cluster_id; ++j)
            {
                if (flags[j])
                    continue;

                double d = DistanceMedian(kmer_vectors[i], kmer_vectors[j]);
                if (d < minimum)
                {
                    c1 = i;
                    c2 = j;
                    minimum = d;
                }
            }
        }

        kmer_vectors[max_cluster_id].swap(kmer_vectors[c1]);
        kmer_vectors[max_cluster_id].insert(kmer_vectors[max_cluster_id].end(), kmer_vectors[c2].begin(), kmer_vectors[c2].end());
        childs[max_cluster_id].push_back(c1);
        childs[max_cluster_id].push_back(c2);
        min_distance[max_cluster_id] = minimum;
        flags[c1] = flags[c2] = 1;

        printf("%d, %d, %d, %4.2f\n", c1, c2, max_cluster_id, minimum);
        fflush(NULL);
        //cout << c1 << " " << c2 << " " << max_cluster_id << " " << minimum << endl;

        ++max_cluster_id;
    }

    return 0;
}

