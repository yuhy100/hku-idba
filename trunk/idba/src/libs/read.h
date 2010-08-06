#ifndef __READ_H_

#define __READ_H_

#include "globals.h"

struct Sequence;

struct Read
{
    unsigned char compressed[22];
    short length;

    unsigned GetNucleotide(unsigned index)
    {
        return (compressed[index>>2] >> ((index&3) << 1)) & 3;
    }
    void SetContent(Sequence *seq);
};

#endif
