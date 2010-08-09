#include "globals.h"
#include "KmerVector.h"
#include "Reader.h"
#include "Utils.h"

#include <omp.h>

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <vector>

using namespace std;

bool is_ignore_first = true;
int num_group = 0;
vector<vector<KmerVector> > groups;

omp_lock_t dist_lock;

double Distance(KmerVector *kv1, KmerVector *kv2)
{
    double sum = 0;
    for (int i = 0; i < kv1->size; ++i)
        sum += abs(kv1->ranks[i] - kv2->ranks[i]);
    return sum;
}

double DistanceMin(vector<KmerVector> &kv1, vector<KmerVector> &kv2)
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
    return tmp.front();
}

double DistanceMax(vector<KmerVector> &kv1, vector<KmerVector> &kv2)
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
    return tmp.back();
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

    double sum = 0;
    for (unsigned i = 0; i < tmp.size(); ++i)
    {
        sum += tmp[i];
    }
    return sum / tmp.size();
    
    sort(tmp.begin(), tmp.end());
    return tmp[tmp.size()/2];
}

int main(int argc, char *argv[])
{
    AddParameter("ignore", &is_ignore_first, SIMPLE);
    ProcessParameters(argc, argv);

    omp_init_lock(&dist_lock);

    num_group = argc - 1;

    groups.resize(num_group);
    for (int i = 0; i < num_group; ++i)
    {
        FastAReader reader(argv[i+1]);
        Sequence seq;
        string comment;

        if (is_ignore_first)
            reader.Read(seq, comment);

        while (reader.Read(seq, comment))
        {
            seq.Encode();
            KmerVector kv;
            kv.Compute(&seq);
            groups[i].push_back(kv);
        }
    }

    for (int i = 0; i < num_group; ++i)
    {
        for (int j = 0; j < num_group; ++j)
        {
            cout << DistanceMin(groups[i], groups[j]) << " ";
        }
        cout << endl;
    }
    cout << endl;

    for (int i = 0; i < num_group; ++i)
    {
        for (int j = 0; j < num_group; ++j)
        {
            cout << DistanceMax(groups[i], groups[j]) << " ";
        }
        cout << endl;
    }
    cout << endl;

    for (int i = 0; i < num_group; ++i)
    {
        for (int j = 0; j < num_group; ++j)
        {
            cout << DistanceMedian(groups[i], groups[j]) << " ";
        }
        cout << endl;
    }
    cout << endl;

    omp_destroy_lock(&dist_lock);

    return 0;
}
