#ifndef __UTILS_H_

#define __UTILS_H_

#include "globals.h"

#include <cstdio>

struct Sequence;
struct HashGraph;
struct ContigGraph;

enum ParameterType { SIMPLE, INTEGER, FLOAT, STRING, };

void AddParameter(const char *name, void *pointer, ParameterType type);
void ProcessParameters(int &argc, char *argv[]);

FILE *OpenFile(const char *filename, const char *mode);

bool ReadFastq(std::FILE *fp, Sequence *seq,
               char *comment = NULL, char *quality = NULL);
bool ReadFasta(std::FILE *fp, Sequence *seq, char *comment = NULL);
void WriteFasta(std::FILE *fp, Sequence *seq, const char *comment = NULL);

uint64 CountSequences(std::FILE *fp, bool isCheck = false);
uint64 ReadHashGraph(std::FILE *fp,
                   HashGraph *hashGraph, unsigned minCount = 0);
uint64 WriteHashGraph(std::FILE *fp,
                      HashGraph *hashGraph, unsigned minCount = 0);

// uint64 ReadContigGraph(std::FILE *fp, ContigGraph *contigGraph);

#endif
