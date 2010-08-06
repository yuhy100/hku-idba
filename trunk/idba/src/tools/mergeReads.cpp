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
char comment2[MaxLine];
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
        fprintf(stderr, "usage: mergeReads read-file1 read-file2 merge-read-file\n");
        exit(1);
    }

    FILE *freadFile1 = OpenFile(argv[1], "rb");
    FILE *freadFile2 = OpenFile(argv[2], "rb");
    FILE *fout = OpenFile(argv[3], "wb");

    Sequence seq;
    Sequence seq2;
    unsigned index = 0;

    while (ReadFasta(freadFile1, &seq, comment))
    {
        if (!ReadFasta(freadFile2, &seq2, comment2))
            break;

        if (length == 0)
        {
            fprintf(fout, ">read%d/1\n", index);
            WriteFasta(fout, &seq);
            fprintf(fout, ">read%d/2\n", index);
            WriteFasta(fout, &seq2);
            ++index;
        }

    }


    fclose(freadFile1);
    fclose(freadFile2);
    fclose(fout);

    return 0;
}
