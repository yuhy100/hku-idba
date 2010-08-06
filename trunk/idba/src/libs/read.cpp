#include "globals.h"
#include "log.h"
#include "sequence.h"
#include "read.h"

#include <cstdio>
#include <algorithm>
#include <iostream>

using namespace std;

void Read::SetContent(Sequence *seq)
{
    if (!seq->IsCodon())
    {
        LogError("CompactSequence::SetSequence: not in codon format\n");
        exit(1);
    }

    for (int i = 0; i < seq->length; ++i)
    {
        compressed[i>>2] |= seq->bases[i] << ((i&3) << 1);
    }
    length = seq->length;
}

