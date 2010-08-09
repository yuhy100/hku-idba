#include "globals.h"
#include "Log.h"
#include "Sequence.h"
#include "Utils.h"
#include "Reader.h"
#include "Writer.h"

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

    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: normReads fa-file norm-fa-file\n");
        fprintf(stderr, "       [--length l] [--mate]\n");
        exit(1);
    }

    FastAReader reader(argv[1]);
    FastAWriter writer(argv[2]);

    Sequence seq;
    Sequence seq2;
    unsigned index = 0;
    string comment;
    string comment2;

    if (mate)
    {
        while (reader.Read(seq, comment))
        {
            if (!reader.Read(seq2, comment2))
                break;

            if (length == 0)
            {
                writer.Write(seq, comment);
                writer.Write(seq2, comment2);
                index += 2;
            }
            else
            {
                if (seq.Size() >= length && seq2.Size() >= length)
                {
                    seq.Resize(length);
                    seq2.Resize(length);
                    writer.Write(seq, comment);
                    writer.Write(seq2, comment2);
                    index += 2;
                }
            }
        }
    }
    else
    {
        while (reader.Read(seq, comment))
        {
            if (length == 0)
            {
                writer.Write(seq, comment);
                ++index;
            }
            else
            {
                if (seq.Size() >= length)
                {
                    seq.Resize(length);
                    writer.Write(seq, comment);
                    ++index;
                }
            }
        }
    }

//    fclose(freadFile);
//    fclose(ftableFile);

    return 0;
}
