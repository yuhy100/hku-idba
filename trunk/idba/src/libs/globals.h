#ifndef __GLOBALS_H_

#define __GLOBALS_H_

#define DEBUG           1

#define KMERSIZE        2

typedef unsigned long long uint64;

extern int kmerLength;

const uint64 MaxLine = 5000000;
const uint64 MaxSeqLength = 300000000;
const uint64 MaxHashTable = 9999997;
const uint64 MaxReadLength = 80;
const uint64 HashNodeBufferSize = (1 << 15);
const uint64 DefaultContigsNum = (1 << 15);
const uint64 MaxDistance = 1000000;
const uint64 TimeLimit = 10000;

#endif
