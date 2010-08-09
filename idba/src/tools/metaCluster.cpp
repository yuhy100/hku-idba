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

const string Version = "2.0.1";

const int MaxReads = 400000;
const int Times = 40;
const double minDist = 1e100;
const double p = 3.0;

struct Node
{
    Node *parent;
    vector<Node *> childs;
    vector<string> refs;
    string taxon;
    bool isExist;
    KmerVector center;
    vector<KmerVector> points;
};


int aux[MaxLine];
double dist[MaxLine];
int matrix[100][100];


bool isFasta = false;
bool isClassified = false;
bool isValidate = false;
int mapping[MaxLine];
int flags[MaxLine];
double sd[MaxLine];
char line[MaxLine];
int classes = 2;

KmerVector kmerVectors[MaxReads];
int type[MaxReads];
int group[MaxReads];
KmerVector center[MaxReads];
double distSum[MaxReads];
double distSumSquare[MaxReads];
double distMean[MaxReads];
double distSD[MaxReads];
double distNum[MaxReads];
double tmp[MaxReads];
string comments[MaxReads];
Sequence reads[MaxReads];
int best[MaxReads];
int isOutlier[MaxReads];
int total;

string taxon;
int refTaxonID[MaxReads];
int taxonID[MaxReads] = {0};
int flagsType[MaxReads] = {0};
int flagsGroup[MaxReads] = {0};
int assign[100];
int predictID[100];

const int MaxNodes = 1000000;

char buf[MaxLine];
char name[MaxLine];
char attribute[MaxLine];
Node nodes[MaxNodes];
map<string, vector<int> > taxonDict;
string taxonNames[] =
{
    "phylum", "class", "order", "family", "genus", "species",
};

vector<string> refList;
vector<string> allRefs;
vector<KmerVector> allKmerVectors;
map<string, int> refRank;

double Distance(KmerVector *kv1, KmerVector *kv2)
{
    double sum = 0;
    for (int i = 0; i < kv1->size; ++i)
        sum += abs(kv1->ranks[i] - kv2->ranks[i]);
    return sum;
}

void Vote()
{
    for (int i = 0; i < classes; ++i)
    {
        for (int j = 0; j < classes; ++j)
            matrix[i][j] = 0;
    }
    for (int i = 0; i < total; ++i)
    {
        if (!flagsType[type[i]] && !isOutlier[i] && !flagsGroup[group[i]])
        {
            ++matrix[type[i]][group[i]];
        }
    }
}

void GetLabel(int &t, int &l)
{
    t = -1;
    l = -1;
    for (int i = 0; i < classes; ++i)
    {
        for (int j = 0; j < classes; ++j)
        {
            if (!flagsType[i] && !flagsGroup[j] && (t == -1 || matrix[i][j] > matrix[t][l]))
            {
                t = i;
                l = j;
            }
        }
    }

}

double DistributePoints()
{
    double power = 0;
    for (int i = 0; i < classes; ++i)
    {
        distSum[i] = 0;
        distSumSquare[i] = 0;
        distNum[i] = 0;
    }

    for (int i = 0; i < total; ++i)
    {
        if (isOutlier[i])
            continue;

        double minimum = 1e100;
        int t = -1;
        for (int j = 0; j < classes; ++j)
        {
            double d = Distance(&kmerVectors[i], &center[j]);
            if (d < minimum)
            {
                minimum = d;
                t = j;
            }
        }

        if (minimum < minDist)
        {
            type[i] = t;
            power += minimum;
            distSum[t] += minimum;
            distSumSquare[t] += minimum * minimum;
            distNum[t]++;
        }
        else
            type[i] = -1;
    }

    for (int i = 0; i < classes; ++i)
    {
        distMean[i] = distSum[i] / distNum[i];
        distSD[i] = sqrt((distSumSquare[i] - distSum[i]*distMean[i]) / (distNum[i] - 1));
    }

    int count = 0;
    for (int i = 0; i < total; ++i)
    {
        if (isOutlier[i] != true)
            ++count;
    }

    return power / count;
}

bool CalCenter()
{
    bool isUpdate = false;
    for (int j = 0; j < classes; ++j)
    {
        KmerVector sum;
        int count = 0;
        for (int i = 0; i < total; ++i)
        {
            if (!isOutlier[i] && type[i] == j)
            {
                for (int l = 0; l < sum.size; ++l)
                    sum.values[l] += kmerVectors[i].values[l];
                ++count;
            }
        }

        for (int l = 0; l < sum.size; ++l)
            sum.values[l] /= count;

        for (int l = 0; l < center[j].size; ++l)
        {
            if (fabs(center[j].values[l] - sum.values[l]) > 1e-6)
                isUpdate = true;
            center[j].values[l] = sum.values[l];
        }

        center[j].ComputeRank();
    }

    return isUpdate;
}

double Kmeans()
{
    // Randomly select the cluster center.
    for (int i = 0; i < total; ++i)
    {
        aux[i] = i;
        isOutlier[i] = false;
        type[i] = -1;
    }

    for (int i = 0; i < classes; ++i)
    {
        int j = i + rand()%(total-i);
        swap(aux[i], aux[j]);
        center[i] = kmerVectors[aux[i]];
    }

    double power = 0;
    for (int k = 0; k < Times; ++k)
    {
        power = DistributePoints();

        if (CalCenter() == false)
            break;
    }

    power = DistributePoints();

    int count = 0;
    for (int i = 0; i < total; ++i)
    {
        if (isOutlier[i])
            continue;

        double d = Distance(&kmerVectors[i], &center[type[i]]);
        if (d > distMean[type[i]] + 2 * distSD[type[i]])
        {
            isOutlier[i] = true;
            ++count;
        }
    }

    power = 0;
    for (int k = 0; k < Times; ++k)
    {
        power = DistributePoints();

        if (CalCenter() == false)
            break;
    }

    return power;
}

void Normalize()
{
    int size = kmerVectors[0].size;
    for (int i = 0; i < size; ++i)
    {
        double sum = 0;
        double sumSquare = 0;
        for (int j = 0; j < total; ++j)
        {
            sum += kmerVectors[j].values[i];
            sumSquare += kmerVectors[j].values[i] * kmerVectors[j].values[i];
        }
        double mean = sum / total;
        double sd = sqrt((sumSquare - sum*mean) / (total - 1));
        for (int j = 0; j < total; ++j)
        {
            kmerVectors[j].values[i] -= mean;
            kmerVectors[j].values[i] /= sd;
        }
    }
}

int Classify(int t)
{
    vector<int> &targets = taxonDict[taxon];
    refRank.clear();

    for (int i = 0; i < total; ++i)
    {
        if (!isOutlier[i] && type[i] == t)
        {
            int min_index = 0;
            double minimum = 1e100;
            for (unsigned j = 0; j < allRefs.size(); ++j)
            {
                if (!flags[j])
                    continue;

                double d = Distance(&allKmerVectors[j], &kmerVectors[i]);
                if (d < minimum)
                {
                    min_index = j;
                    minimum = d;
                }
            }

            refRank[allRefs[min_index]] += 1;
        }
    }

    int r = 0;
    int maximum = 0;
    for (int i = 0; i < targets.size(); ++i)
    {
        vector<string> &refs = nodes[targets[i]].refs;
        int d = 0;

        for (int j = 0; j < refs.size(); ++j)
            d += refRank[refs[j]];

        if (d > maximum)
        {
            maximum = d;
            r = i;
        }
    }

    return targets[r];
}


void ReadTaxonTreeDB(char *db_file)
{
    FILE *fp = OpenFile(db_file, "rb");

    int total;
    fscanf(fp, "%d %s", &total, buf);

    allRefs.resize(total);
    allKmerVectors.resize(total);
    for (int i = 0; i < total; ++i)
    {
        fscanf(fp, "%s", buf);
        allRefs[i] = buf;
        for (int j = 0; j < allKmerVectors[i].size; ++j)
        {
            fscanf(fp, "%lf", &allKmerVectors[i].values[j]);
            allKmerVectors[i].ComputeRank();
        }
    }


    for (int i = 0; i < 6; ++i)
    {
        fscanf(fp, "%s", buf);
        vector<int> &targets = taxonDict[buf];

        int num_group;
        fscanf(fp, "%d", &num_group);
        targets.resize(num_group);
        for (int j = 0; j < num_group; ++j)
        {
            int id, size;
            fscanf(fp, "%d %d", &id, &size);
            targets[j] = id;
            nodes[id].refs.resize(size);
            for (int k = 0; k < size; ++k)
            {
                fscanf(fp, "%s", buf);
                nodes[id].refs[k] = buf;
            }
        }
    }
}

void Usage()
{
    printf("metaCluster %s\n", Version.c_str());
    printf("Usage: metacluster read-file [taxon_tree.db] [--fasta] [--classes c] [--classify]\n");
}

int main(int argc, char *argv[])
{
    srand((unsigned)time(NULL));

    AddParameter("classes", &classes, INTEGER);
    AddParameter("classify", &isClassified, SIMPLE);
    AddParameter("fasta", &isFasta, SIMPLE);
    AddParameter("validate", &isValidate, SIMPLE);
    ProcessParameters(argc, argv);

    if (argc < 2)
    {
        //fprintf(stderr, "usage: metaCluster read-file\n");
        Usage();
        exit(1);
    }

    // Get the reads and compute the kmer vectors
    Sequence seq;
    unsigned index = 0;
    FastAReader reader(argv[1]);
    string comment;
    while (reader.Read(seq, comment))
    {
            if (seq.IsChar())
            {
                comments[index] = comment;
                reads[index] = seq;
                seq.Encode();
                kmerVectors[index].Compute(&seq);

                ++index;

                strcpy(line, comment.c_str());
                if (isValidate)
                {
                    for (int i = 0; line[i]; ++i)
                    {
                        if (line[i] == '_')
                        {
                            group[index] = atoi(line + i+1);
                            break;
                        }
                    }
                }
            }
    }

    total = index;

//     Normalize();


    double optimal = 1e100;
    for (int round = 0; round < 10; ++round)
    {
        double power = Kmeans();

        if (power < optimal)
        {
            optimal = power;
            copy(type, type + total, best);
        }
    }

    copy(best, best + total, type);

    int correct = 0;
    fill_n(flagsType, classes, 0);
    fill_n(flagsGroup, classes, 0);

    if (isValidate)
    {
        for (int i = 0; i < classes; ++i)
        {
            Vote();
            int t = 0, l = 0;
            GetLabel(t, l);
            assign[t] = l;
            flagsType[t] = 1;
            flagsGroup[l] = 1;

            int count = 0;
            int corr = 0;
            for (int j = 0; j < total; ++j)
            {
                if (!isOutlier[j] && type[j] == t)
                {
                    ++count;
                    if (group[j] == l)
                        ++corr;
                }
            }

            correct += corr;
                    
            printf("%d/%d type %d group %d %.2f%%\n", corr, count, t, l, corr*100.0/count);
        }
    }

    if (isClassified)
    {
        //BuildTaxonTree(argc, argv);
        ReadTaxonTreeDB(argv[2]);

        //cout << "read tree finished" << endl;
        
        taxon = argv[3];
        vector<int> &targets = taxonDict[argv[3]];

        //cout << argv[4] << endl;
        //cout << targets.size() << endl;
        printf("taxon level %s, total %d %s in database.\n", argv[3], (int)targets.size(), argv[3]);

        fill_n(flags, MaxLine, 1);
        if (isValidate)
        {
            FILE *frefs = OpenFile(argv[4], "rb");
            while (fscanf(frefs, "%s", buf) != EOF)
                refList.push_back(buf);
            for (int i = 0; i < refList.size(); ++i)
            {
                refTaxonID[i] = 0;
                for (int j = 0; j < targets.size(); ++j)
                {
                    vector<string> &refs = nodes[targets[j]].refs;
                    if (find(refs.begin(), refs.end(), refList[i]) != refs.end())
                    {
                        refTaxonID[i] = targets[j];
                        break;
                    }
                }
            }
            for (unsigned i = 0; i < allRefs.size(); ++i)
            {
                if (find(refList.begin(), refList.end(), allRefs[i]) != refList.end())
                    flags[i] = false;
            }
        }

        for (int t = 0; t < classes; ++t)
        {
            int m = assign[t];
            int id = Classify(t);
            predictID[t] = id;

            int count = 0;
            for (int i = 0; i < total; ++i)
            {
                if (!isOutlier[i] && type[i] == t)
                {
                    taxonID[i] = id;
                    ++count;
                }
            }

            printf("group %d, %d fragments ", t, count);
            printf("assigned id: %d ", id);

            if (isValidate)
                printf(", original taxon id: %d\n", refTaxonID[m]);
            else
                printf("\n");
        }
    }
    else
    {
        for (int t = 0; t < classes; ++t)
        {
            int count = 0;
            for (int i = 0; i < total; ++i)
            {
                if (!isOutlier[i] && type[i] == t)
                    ++count;
            }

            printf("group %d, %d fragments\n", t, count);
        }
    }

    if (isValidate)
    {
        int count = 0;
        int classCount = 0;
        for (int i = 0; i < total; ++i)
        {
            if (isOutlier[i])
                ++count;
            else if (taxonID[i] == refTaxonID[group[i]])
                ++classCount;
        }

        //printf("%.2f%% ", correct*100.0/(total - count));
        printf("clustering accuracy: %d/%d %.2f%%\n", correct, total - count, correct*100.0/(total - count));
        printf("classification accuracy: %.2f%%\n", classCount*100.0/(total - count));
    }

    FILE *foutput[100];
    for (int i = 0; i < classes; ++i)
    {
        sprintf(buf, "%s-out%03d.fa", argv[1], i);
        foutput[i] = OpenFile(buf, "wb");
        if (isClassified)
            fprintf(foutput[i], ">taxon id %d\n", predictID[i]);
        fflush(NULL);
    }

    for (int i = 0; i < total; ++i)
    {
        if (!isOutlier[i])
        {
            FILE *fp = foutput[type[i]];
            fprintf(fp, ">%s\n", comments[i].c_str());

            if (isFasta)
            {
                fprintf(fp, "%s\n", reads[i].ToCString());
            }
        }
    }

    return 0;
}

