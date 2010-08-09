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
    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: fq2fa fq-file fa-file\n");
        exit(1);
    }

    FastQReader reader(argv[1]);
    FastAWriter writer(argv[2]);
    
    Sequence seq;
    string ss;
    int index = 0;
    while (reader.Read(seq, ss))
    {
        writer.Write(seq, ss);
        ++index;
    }
    printf("total %d reads\n", index);

    return 0;
}
