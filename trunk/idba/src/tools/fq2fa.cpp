#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "utils.h"
#include "hashGraph.h"
#include "read.h"

#include <cstdio>
#include <algorithm>
#include <cstring>

using namespace std;

char line[MaxLine];
char comment[MaxLine];
FILE *fp[MaxLine];
int length = 0;
bool mate = 0;

int main(int argc, char *argv[])
{
    ProcessParameters(argc, argv);

    HashGraph hashGraph(MaxHashTable);
    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: fq2fa fq-file fa-file\n");
        exit(1);
    }

    FILE *freadFile = OpenFile(argv[1], "rb");
    FILE *ftableFile = OpenFile(argv[2], "wb");

    Sequence seq;
    Sequence seq2;
    unsigned index = 0;

    while (ReadFastq(freadFile, &seq, comment))
    {
        WriteFasta(ftableFile, &seq, comment);
        ++index;
    }

    printf("total %d reads\n", index);

    fclose(freadFile);
    fclose(ftableFile);

    return 0;
}
