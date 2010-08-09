#include "globals.h"
#include "Log.h"
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

const int MaxReads = 400000;
const int MaxNodes = 1000000;

char line[MaxLine];
char buf[MaxLine];
char name[MaxLine];
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

void dfs(Node *node)
{
    if (node->childs.size() > 0)
    {
        int count = 0;
        node->center.Clear();
        for (unsigned i = 0; i < node->childs.size(); ++i)
        {
            dfs(node->childs[i]);
            if (node->childs[i]->refs.size() > 0)
            {
                for (unsigned j = 0; j < node->childs[i]->refs.size(); ++j)
                {
                    node->refs.push_back(node->childs[i]->refs[j]);
                    node->points.push_back(node->childs[i]->points[j]);
                }
                node->center += node->childs[i]->center;
                ++count;
            }
        }

        if (count != 0)
            node->center /= count;
    }
}

void BuildTaxonTree(char *nodes_file, char *gff_list)
{
    FILE *fnodes = OpenFile(nodes_file, "rb");
    FILE *fgff = OpenFile(gff_list, "rb");

    for (int i = 0; i < MaxNodes; ++i)
        nodes[i].isExist = false;

    while (fgets(line, MaxLine, fnodes) != NULL)
    {
        int id, pid;
        sscanf(line, "%d %s %d %s %s", &id, buf, &pid, buf, name);
        nodes[id].parent = &nodes[pid];
        nodes[id].taxon = name;
        nodes[id].isExist = true;

        if (id != pid)
            nodes[pid].childs.push_back(&nodes[id]);
    }

    fclose(fnodes);

    while (fscanf(fgff, "%s", name) != EOF)
    {
        FILE *fp = OpenFile(name, "rb");

        while (fgets(line, MaxLine, fp) != NULL)
        {
            if (line[0] == '#')
                continue;
            else
            {
                if (strstr(line, "plasmid") == NULL)
                {
                    char *p = strstr(line, "taxon");

                    if (p == NULL)
                    {
                        printf("%s\n", name);
                        fprintf(stderr, "%s\n", p);
                        fprintf(stderr, "gff format error\n");
                        exit(1);
                    }

                    int id = atoi(p + 6);

                    if (nodes[id].isExist == false)
                    {
                        fprintf(stderr, "id not exist\n");
                        exit(1);
                    }

                    strcpy(name + strlen(name)-3, "fna");

                    bool flag = false;
                    for (unsigned i = 0; i < refList.size(); ++i)
                    {
                        if (name == refList[i])
                        {
                            flag = true;
                            break;
                        }
                    }

                    if (flag)
                        break;

                    nodes[id].refs.push_back(name);

                    FILE *fref = OpenFile(name, "rb");
                    Sequence seq;
                    //ReadFasta(fref, &seq);

                    FastAReader reader(name);
                    string comment;
                    reader.Read(seq, comment);
                    seq.Encode();
                    nodes[id].center.Compute(&seq);
                    nodes[id].points.push_back(nodes[id].center);
                    fclose(fref);

                    allRefs.push_back(name);
                    allKmerVectors.push_back(nodes[id].center);
                }

                break;
            }
        }

        fclose(fp);
    }

    fclose(fgff);

    for (int i = 0; i < MaxNodes; ++i)
    {
        if (nodes[i].isExist)
        {
            Node *node = nodes[i].parent;
            while (node != node->parent)
            {
                bool found = false;
                for (int j = 0; j < 6; ++j)
                {
                    if (node->taxon == taxonNames[j])
                    {
                        found = true;
                        break;
                    }
                }

                if (found)
                    break;

                node = node->parent;
            }

            nodes[i].parent = node;
        }
    }

    dfs(&nodes[1]);

    for (int i = 0; i < MaxNodes; ++i)
    {
        if (nodes[i].isExist && nodes[i].refs.size() > 0)
        {
            taxonDict[nodes[i].taxon].push_back(i);
        }
    }

//    taxon = argv[5];
//    vector<int> &targets = taxonDict[argv[5]];
//
//    for (int i = 0; i < refList.size(); ++i)
//    {
//        strcpy(name, refList[i].c_str());
//        strcpy(name + strlen(name) - 3, "gff");
//
//        FILE *fp = OpenFile(name, "rb");
//
//        while (fgets(line, MaxLine, fp) != NULL)
//        {
//            if (line[0] == '#')
//                continue;
//            else
//            {
//                if (strstr(line, "plasmid") == NULL)
//                {
//                    char *p = strstr(line, "taxon");
//
//                    if (p == NULL)
//                    {
//                        printf("%s\n", name);
//                        fprintf(stderr, "%s\n", p);
//                        fprintf(stderr, "gff format error\n");
//                        exit(1);
//                    }
//
//                    int id = atoi(p + 6);
//
//                    if (nodes[id].isExist == false)
//                    {
//                        fprintf(stderr, "id not exist\n");
//                        exit(1);
//                    }
//
//                    Node *parent = nodes[id].parent;
//                    while (parent->taxon != taxon && parent->parent != parent)
//                        parent = parent->parent;
//                    refTaxonID[i] = parent - nodes;
//                }
//
//                break;
//            }
//        }
//
//        fclose(fp);
//    }
}

bool isClassified = false;

int main(int argc, char *argv[])
{
    ProcessParameters(argc, argv);

    if (argc < 3)
    {
        fprintf(stderr, "usage: cluster read-file ref-list\n");
        exit(1);
    }

    BuildTaxonTree(argv[1], argv[2]);

    printf("%d references\n", (int)allRefs.size());
    for (unsigned i = 0; i < allRefs.size(); ++i)
    {
        printf("%s\n", allRefs[i].c_str());
        for (int j = 0; j < allKmerVectors[i].size; ++j)
        {
            printf("%.6f ", allKmerVectors[i].values[j]);
        }
        printf("\n");
    }
    printf("\n");

    for (int i = 0; i < 6; ++i)
    {
        vector<int> &targets = taxonDict[taxonNames[i]];
        printf("%s %d\n", taxonNames[i].c_str(), (int)targets.size());
        for (unsigned j = 0; j < targets.size(); ++j)
        {
            vector<string> &refs = nodes[targets[j]].refs;
            printf("%d %d\n", targets[j], (int)refs.size());
            for (unsigned k = 0; k < refs.size(); ++k)
                printf("%s\n", refs[k].c_str());
        }
        printf("\n");
    }

    return 0;
}

