#include "globals.h"
#include "Log.h"
#include "Sequence.h"
#include "Utils.h"
#include "Reader.h"
#include "Writer.h"

#include <cctype>
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
//char genome[1000000000];
int errorNums[1000] = {0};
bool separate = false;


int main(int argc, char *argv[])
{
    AddParameter("depth", &depth, INTEGER);
    AddParameter("errorRate", &errorRate, FLOAT);
    AddParameter("readLength", &readLength, INTEGER);
//     AddParameter("readNumber", &readNumber, INTEGER);
    AddParameter("mate", &mate, SIMPLE);
    AddParameter("insertLength", &insertLength, INTEGER);
    AddParameter("separate", &separate, SIMPLE);

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


    FastAReader reader(argv[1]);
    FastAWriter writer1(argv[2]);

    //int genomeLength = 0;
    Sequence ref;
    string comment;
    //while (ReadFasta(fref, &ref))
    Sequence genome;
    while (reader.Read(ref, comment))
    {
        //copy(ref.ToCString(), ref.ToCString() + ref.Size(), genome + genomeLength);
        if (genome.Size() != 0)
            genome += 'N';
        genome += ref;
//        genomeLength += ref.Size();
//        for (int i = 0; i < insertLength; ++i)
//            genome[genomeLength++] = 'N';
    }
    //ref.SetContent(genome, genomeLength);
    //genome[genomeLength] = '\0';
    ref = genome;

    for (unsigned i = 0; i < ref.Size(); ++i)
    {
        ref[i] = toupper(ref[i]);
    }

    if (readNumber == 0)
        readNumber = (depth * ref.Size()) / readLength;
    else
        depth = int(readLength*readNumber/ref.Size());

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
                offset = rand() % (ref.Size() - readLength + 1);
                ref.GetSubSequence(seq, offset, readLength);
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
                    int c = seq[j];
                    while (c == seq[j])
                    {
                        c = rand()/93 % 4;
                    }
                    seq[j] = c;
                }
            }
            seq.Decode();

            ++errorNums[error];

            if (rand() < RAND_MAX/2)
                seq.ReverseComplement();

//            fprintf(freadFile, ">read%d\n", i);
//            WriteFasta(freadFile, &seq);
            writer1.WriteFormat(seq, "read%d", i);
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
                offset = rand() % (ref.Size() + 1 - insertLength);
                ref.GetSubSequence(seq1, offset, readLength);
                ref.GetSubSequence(seq2, offset + insertLength - readLength, readLength);
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
                    int c = seq1[j];
                    while (c == seq1[j])
                    {
                        c = rand()/93 % 4;
                    }
                    seq1[j] = c;
                }

                if (rand()*1.0 < errorRate*RAND_MAX)
                {
                    ++error2;
                    int c = seq2[j];
                    while (c == seq2[j])
                    {
                        c = rand()/93 % 4;
                    }
                    seq2[j] = c;
                }
            }

            ++errorNums[error1];
            ++errorNums[error2];

            seq1.Decode();
            seq2.Decode();
            seq2.ReverseComplement();

            if (i%2 == 0)
            {
                writer1.WriteFormat(seq2, "read%d/1", i);
                writer1.WriteFormat(seq1, "read%d/2", i);
            }
            else
            {
                writer1.WriteFormat(seq1, "read%d/1", i);
                writer1.WriteFormat(seq2, "read%d/2", i);
            }
        }
    }


    return 0;
}
