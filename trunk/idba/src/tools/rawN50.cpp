#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "utils.h"
#include "compactSequence.h"

#include <cctype>
#include <cstdio>
#include <algorithm>
#include <cstring>

using namespace std;

const int MaxContigs = 1000000;

Sequence contigs[MaxContigs];
int len[MaxContigs];
char line[MaxLine];
int minContig = 100;
int refLength = 0;

int main(int argc, char *argv[])
{
    AddParameter("minContig", &minContig, INTEGER);
    AddParameter("refLength", &refLength, INTEGER);
    ProcessParameters(argc, argv);



    Sequence seq;
    unsigned index = 0;
    while (ReadFasta(stdin, &seq))
    {
        {
            if (seq.length < minContig)
                continue;

            for (int i = 0; i < seq.length; ++i)
                seq.bases[i] = toupper(seq.bases[i]);

            if (!seq.IsChar())
                continue;

            contigs[index].SetContent(&seq);
            ++index;
        }
    }

    int total = index;
    int sum = 0;
    for (int i = 0; i < total; ++i)
    {
        len[i] = contigs[i].length;
        sum += len[i];
    }

    sort(len, len + total);
    reverse(len, len + total);

    int n50 = 0;
    int s = 0;

    if (refLength > 0)
        sum = refLength;
    int n80 = 0;
    for (int i = 0; i < total; ++i)
    {
        s += len[i];
        if (s > sum*0.5 && n50 == 0)
            n50 = len[i];
        if (s > sum*0.8 && n80 == 0)
            n80 = len[i];
    }

    printf("contigs: %d n50: %d max: %d mean: %d total length: %d n80: %d\n\n",
           total, n50, len[0], sum/total, sum, n80);

    return 0;
}
