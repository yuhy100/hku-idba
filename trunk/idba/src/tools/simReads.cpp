#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "utils.h"

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <iostream>

using namespace std;

bool mate = false;
int insertLength = 250;
char line[MaxLine];
uint64 depth = 100;
double errorRate = 0.01;
int readLength = 35;
uint64 readNumber = 0;
char genome[1000000000];
int errorNums[1000] = {0};

int main(int argc, char *argv[])
{
    AddParameter("depth", &depth, INTEGER);
    AddParameter("errorRate", &errorRate, FLOAT);
    AddParameter("readLength", &readLength, INTEGER);
//     AddParameter("readNumber", &readNumber, INTEGER);
    AddParameter("mate", &mate, SIMPLE);
    AddParameter("insertLength", &insertLength, INTEGER);

    ProcessParameters(argc, argv);

    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: simRead ref-file read-file\n");
        fprintf(stderr, "       [--depth d] [--mate]\n");
        fprintf(stderr, "       [--errorRate e] [--readLength l]\n");
        exit(1);
    }

    FILE *fref = fopen(argv[1], "rb");
    if (fref == NULL)
    {
        fprintf(stderr, "open %s failed\n", argv[1]);
        exit(1);
    }

    FILE *freadFile = fopen(argv[2], "wb");
    if (freadFile == NULL)
    {
        fprintf(stderr, "open %s failed\n", argv[2]);
        exit(1);
    }

    int genomeLength = 0;
    Sequence ref;
    while (ReadFasta(fref, &ref))
    {
        copy(ref.bases, ref.bases + ref.length, genome + genomeLength);
        genomeLength += ref.length;
        genome[genomeLength++] = 'N';
    }
    ref.SetContent(genome, genomeLength);

    if (readNumber == 0)
        readNumber = (depth * ref.length) / readLength;
    else
        depth = int(readLength*readNumber/ref.length);

    readNumber += readNumber % 2;
    LogMessage("depth %d\treadLength %d\terrorRate %.3f\treadNumber %d\n",
        depth, readLength, errorRate, readNumber);
    LogMessage("insertLength %d\n", insertLength);

    if (!mate)
    {
        Sequence seq;
        for (unsigned i = 0; i < readNumber; ++i)
        {
            unsigned offset = 0;
            while (true)
            {
                offset = rand() % (ref.length - readLength + 1);
                ref.GetSubSequence(&seq, offset, readLength);
                if (seq.IsChar())
                    break;
            }

            seq.Encode();

            int error = 0;
            for (int j = 0; j < readLength; ++j)
            {
                if (rand()*1.0 < errorRate*RAND_MAX)
                {
                    ++error;
                    int c = seq.bases[j];
                    while (c == seq.bases[j])
                    {
                        c = rand()/93 % 4;
                    }
                    seq.bases[j] = c;
                }
            }
            seq.Decode();

            ++errorNums[error];

            if (rand() < RAND_MAX/2)
                seq.ReverseComplement();

            fprintf(freadFile, ">read%d\n", i);
            WriteFasta(freadFile, &seq);
        }
    }
    else
    {
        Sequence seq1;
        Sequence seq2;
        for (unsigned i = 0; i < readNumber; i += 2)
        {
            unsigned offset = 0;
            while (true)
            {
                offset = rand() % (ref.length + 1 - insertLength);
                ref.GetSubSequence(&seq1, offset, readLength);
                ref.GetSubSequence(&seq2, offset + insertLength - readLength, readLength);
                if (seq1.IsChar() && seq2.IsChar())
                    break;
            }

            seq1.Encode();
            seq2.Encode();
            int error1 = 0;
            int error2 = 0;
            for (int j = 0; j < readLength; ++j)
            {
                if (rand()*1.0 < errorRate*RAND_MAX)
                {
                    ++error1;
                    int c = seq1.bases[j];
                    while (c == seq1.bases[j])
                    {
                        c = rand()/93 % 4;
                    }
                    seq1.bases[j] = c;
                }

                if (rand()*1.0 < errorRate*RAND_MAX)
                {
                    ++error2;
                    int c = seq2.bases[j];
                    while (c == seq2.bases[j])
                    {
                        c = rand()/93 % 4;
                    }
                    seq2.bases[j] = c;
                }
            }

            ++errorNums[error1];
            ++errorNums[error2];

            seq1.Decode();
            seq2.Decode();
            seq2.ReverseComplement();

            if (i%2 == 0)
            {

                fprintf(freadFile, ">read%d/1\n", i);
                WriteFasta(freadFile, &seq2);

                fprintf(freadFile, ">read%d/2\n", i);
                WriteFasta(freadFile, &seq1);
            }
            else
            {
                fprintf(freadFile, ">read%d/1\n", i);
                WriteFasta(freadFile, &seq1);

                fprintf(freadFile, ">read%d/2\n", i);
                WriteFasta(freadFile, &seq2);
            }
        }
    }

    fclose(fref);
    fclose(freadFile);

//     int sum = 0;
//     for (int i = 0; i <= readLength; ++i)
//     {
//         printf("%d %d %d\n", i, errorNums[i], sum += errorNums[i]);
//     }

    return 0;
}
