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
    AddParameter("length", &length, INTEGER);
    AddParameter("mate", &mate, SIMPLE);

    ProcessParameters(argc, argv);

    HashGraph hashGraph(MaxHashTable);
    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: normReads fa-file norm-fa-file\n");
        fprintf(stderr, "       [--length l] [--mate]\n");
        exit(1);
    }

    FILE *freadFile = OpenFile(argv[1], "rb");
    FILE *ftableFile = OpenFile(argv[2], "wb");

    Sequence seq;
    Sequence seq2;
    unsigned index = 0;

    if (mate)
    {
        while (ReadFasta(freadFile, &seq, comment))
        {
            if (!ReadFasta(freadFile, &seq2, comment2))
                break;

            if (length == 0)
            {
                WriteFasta(ftableFile, &seq, comment);
                WriteFasta(ftableFile, &seq2, comment2);
                index += 2;
            }
            else
            {
                if (seq.length >= length && seq2.length >= length)
                {
                    seq.Resize(length);
                    seq2.Resize(length);
                    WriteFasta(ftableFile, &seq, comment);
                    WriteFasta(ftableFile, &seq2, comment2);
                    index += 2;
                }
            }
        }
    }
    else
    {
        while (ReadFasta(freadFile, &seq, comment))
        {
            if (length == 0)
            {
                WriteFasta(ftableFile, &seq, comment);
                ++index;
            }
            else
            {
                if (seq.length >= length)
                {
                    seq.Resize(length);
                    WriteFasta(ftableFile, &seq, comment);
                    ++index;
                }
            }
        }
    }

    fclose(freadFile);
    fclose(ftableFile);

    return 0;
}
