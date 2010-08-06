#include "globals.h"
#include "sequence.h"
#include "utils.h"
#include "log.h"

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <cstring>
#include <set>
#include <string>
#include <map>

using namespace std;

struct Interval
{
    int from, to;
};

int minContig = 100;

char buf[MaxLine];
char line[MaxLine];
char readName[MaxLine];
char refName[MaxLine];
char comment[MaxLine];
vector<string> names;
vector<Sequence> sequences;
map<string, int> dict;
vector<char *> matchs;
vector<vector<int> > lengths;
vector<set<string> > founds;
set<string> allFound;

int main(int argc, char *argv[])
{
    AddParameter("minContig", &minContig, INTEGER);

    ProcessParameters(argc, argv);

    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: calcCover contig-file ref-file\n");
        fprintf(stderr, "       [--minContig l]\n");
        exit(1);
    }

    FILE *fcontigs = OpenFile(argv[1], "rb");
    FILE *fref = OpenFile(argv[2], "rb");

    Sequence seq;
    while (ReadFasta(fref, &seq, comment))
    {
        sscanf(comment, "%s", refName);
        dict[refName] = sequences.size();
        sequences.push_back(seq);
        char *match = new char[seq.length];
        fill_n(match, seq.length, 0);
        matchs.push_back(match);
        lengths.push_back(vector<int>());
        names.push_back(refName);
        founds.push_back(set<string>());
    }

    set<string> found;
    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);
    while (fgets(line, MaxLine, stdin) != NULL)
    {
        if (line[0] == '-' && line[1] == '-')
            break;

        int from, to;
        sscanf(line, "%s %s %s %s %s %s %d %d",
               buf, readName, buf, buf, buf, refName, &from, &to);

        if (from > to)
            swap(from, to);

        if (to - from >= minContig)
        {
            if (dict.find(refName) == dict.end())
            {
                LogError("Input contain unknown sequence\n");
                exit(1);
            }

            int index = dict[refName];
            founds[index].insert(readName);
            allFound.insert(readName);
            lengths[index].push_back(to - from);
            for (int i = from; i < to; ++i)
                matchs[index][i] = 1;
        }
    }

    for (unsigned i = 0; i < sequences.size(); ++i)
    {
        vector<int> &v = lengths[i];
        int n50 = 0;
        int sum = 0;
        int n80 = 0;
        int maximum = 0;
        sort(v.begin(), v.end());
        reverse(v.begin(), v.end());
        for (unsigned j = 0; j < v.size(); ++j)
        {
            sum += v[j];
            if (v[j] > maximum)
                maximum = v[j];
            if (sum > 0.5 * sequences[i].length && n50 == 0)
                n50 = v[j];
            if (sum > 0.8 * sequences[i].length && n80 == 0)
                n80 = v[j];
        }

        int total = 0;
        int contigs = 0;
        char *match = matchs[i];
        for (int j = 0; j < sequences[i].length; ++j)
        {
            total += match[j];
            if (match[j] && j+1 < sequences[i].length && match[j+1] == 0)
                ++contigs;
        }

        printf("%s\nfound: %lu n50: %d n80: %d max: %d mean: %lu\n",
               names[i].c_str(), founds[i].size(), n50, n80, maximum, sum/v.size());
        printf("coverage: %.2f%% merged contigs: %d mapped length: %d %d\n",
               total*100.0/sequences[i].length, contigs, sum, total);
    }

    int wrong = 0;
    int wrongLength = 0;
    fseek(fcontigs, 0, SEEK_SET);
    while (ReadFasta(fcontigs, &seq, readName))
    {
        if (seq.length >= minContig && allFound.find(readName) == allFound.end())
        {
            ++wrong;
            wrongLength += seq.length;
        }
    }
    printf("wrong: %d wrong length: %d\n", wrong, wrongLength);

    sprintf(line, "%s-e", argv[2]);
    FILE *fcontigsE = OpenFile(line, "wb");

    while (ReadFasta(fcontigs, &seq, comment))
    {
        if (allFound.find(comment) == allFound.end())
            WriteFasta(fcontigsE, &seq, comment);
    }

    return 0;
}
