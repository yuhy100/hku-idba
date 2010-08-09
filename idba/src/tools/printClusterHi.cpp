#include "globals.h"

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <string>
#include <cctype>

using namespace std;

const int MaxCluster = 2048;

char buf[MaxLine];
char line[MaxLine];
vector<string> filenames;

vector<int> childs[MaxCluster];
double min_distance[MaxCluster] = {0};
double value_print[MaxCluster] = {0};

void ReplaceComma(char *s)
{
    for (int i = 0; s[i]; ++i)
    {
        if (s[i] == ',')
            s[i] = ' ';
    }
}

void PrintCluster(int id)
{
    if (id < filenames.size())
    {
        int num_round = filenames.size() - 1;

        string s;
        for (int i = num_round - 1; i >= 0; --i)
        {
            if (value_print[i + filenames.size()])
                sprintf(buf, "%04.2f---", min_distance[i + filenames.size()]);
            else
                sprintf(buf, "----------");
            s += buf;
        }
        s += filenames[id];
        printf("%s\n", s.c_str());
    }
    else
    {
        value_print[id] = true;
        PrintCluster(childs[id][0]);
        PrintCluster(childs[id][1]);
        value_print[id] = false;
    }
}

int main(int argc, char *argv[])
{
    while (fgets(line, MaxLine, stdin) != NULL)
    {
        if (!isdigit(line[0]))
        {
            sscanf(line, "%s", buf);
            filenames.push_back(buf);
        }
        else
        {
            int c1, c2, c3;
            double minimum;
            ReplaceComma(line);
            sscanf(line, "%d %d %d %lf", &c1, &c2, &c3, &minimum);
            childs[c3].push_back(c1);
            childs[c3].push_back(c2);
            min_distance[c3] = minimum;
        }
    }

    PrintCluster(filenames.size() + filenames.size() - 2);

    return 0;
}

