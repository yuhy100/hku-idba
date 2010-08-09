#include "globals.h"
//#include "Log.h"
#include "Sequence.h"
#include "Read.h"

#include <cstdio>
#include <algorithm>
#include <iostream>

using namespace std;

void Read::SetContent(const Sequence &seq)
{
    for (int i = 0; i < seq.Size(); ++i)
        SetNucleotide(i, seq[i]);
    length = seq.Size();
    status = 0;
}

