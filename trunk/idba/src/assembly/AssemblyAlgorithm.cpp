#include "globals.h"
#include "AssemblyAlgorithm.h"
#include "Reader.h"
#include "Writer.h"
#include "Sequence.h"

#include <iostream>
#include <cstdio>
#include <algorithm>
#include <string>

using namespace std;


int64 AssemblyAlgorithm::LoadReads(const string &readfile)
{
    FastAReader reader(readfile);
    Sequence seq;

    int64 numReads = reader.NumReads();
    reads.resize(numReads);

    int64 index = 0;
    string comment;
    while (reader.Read(seq, comment))
    {
        seq.Trim(trim);
        seq.TrimError();
        if (seq.IsChar() && seq.Size() >= kmerLength)
        {
            seq.Encode();

            if (seq.Size() < (int)kmerLength)
                continue;

            if (seq.Size() <= (int)Read::MaxReadLength)
                reads[index++] = seq;
            else
                longReads.push_back(seq);
        }
        else
            reads[index++].Resize(0);
    }

    return numReads;
}
