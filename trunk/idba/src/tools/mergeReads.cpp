#include "globals.h"
#include "Log.h"
#include "Sequence.h"
#include "Reader.h"
#include "Writer.h"

#include <cstdio>
#include <algorithm>
#include <cstring>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 4
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: mergeReads read-file1 read-file2 merge-read-file\n");
        exit(1);
    }

    FastAReader reader1(argv[1]);
    FastAReader reader2(argv[2]);
    FastAWriter writer(argv[3]);

    Sequence seq1;
    Sequence seq2;
    string comment1;
    string comment2;

    int64 index = 0;
    while (reader1.Read(seq1, comment1) && reader2.Read(seq2, comment2))
    {
        writer.Write(seq1, comment1);
        writer.Write(seq2, comment2);
        ++index;
    }

    cout << index << " pairs." << endl;

    return 0;
}
