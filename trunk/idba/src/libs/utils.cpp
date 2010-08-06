#include "globals.h"
#include "log.h"
#include "utils.h"
#include "sequence.h"
#include "hashGraph.h"
#include "contigGraph.h"

#include <cstdio>
#include <cstring>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <map>
#include <string>

using namespace std;

static const int FastaMaxLine = 80;
static char line[MaxLine];
static char comment[MaxLine];
static char buf[MaxSeqLength];

static map<string, void *> parameters;
static map<string, ParameterType> types;

void AddParameter(const char *name, void *pointer, ParameterType type)
{
    parameters[name] = pointer;
    types[name] = type;
}

void ProcessParameters(int &argc, char *argv[])
{
    int k = 0;
    for (int i = 0; i < argc; ++i)
    {
        if (argv[i][0] == '-' && argv[i][1] == '-')
        {
            char *name = argv[i] + 2;
            if (parameters.find(name) != parameters.end())
            {
                switch (types[name])
                {
                    case SIMPLE:
                        *(bool *)parameters[name] = true;
                        break;

                    case INTEGER:
                        *(int *)parameters[name] = atoi(argv[++i]);
                        break;

                    case FLOAT:
                        *(double *)parameters[name] = atof(argv[++i]);
                        break;

                    case STRING:
                        strcpy((char*)parameters[name], argv[++i]);
                }
            }
        }
        else
        {
            argv[k++] = argv[i];
        }
    }
    argv[k] = NULL;
    argc = k;
}

FILE *OpenFile(const char *filename, const char *mode)
{
    FILE *fp = fopen(filename, mode);

    if (fp == NULL)
    {
        LogError("open %s failed\n", filename);
        exit(1);
    }

    return fp;
}

bool ReadFastq(FILE *fp, Sequence *seq, char *comment, char *quality)
{
    if (fgets(line, MaxLine, fp) == NULL)
        return false;

    if (line[0] != '@')
    {
        LogError("ReadFastq: comment is not in correct format\n");
        exit(1);
    }

    line[strlen(line)-1] = '\0';
    if (comment != NULL)
        strcpy(comment, line + 1);

    int offset = 0;
    while (fgets(line, MaxLine, fp) != NULL)
    {
        if (line[0] == '+')
        {
            fgets(line, MaxLine, fp);
            line[strlen(line)-1] = '\0';
            if (quality != NULL)
                strcpy(quality, line);
            break;
        }
        else
        {
            int len = strlen(line);
            while (isspace(line[len-1]))
                --len;
            copy(line, line + len, buf + offset);
            offset += len;
        }
    }

    seq->SetContent(buf, offset);

    return offset != 0;
}

bool ReadFasta(FILE *fp, Sequence *seq, char *comment)
{
    if (fgets(line, MaxLine, fp) == NULL)
        return false;

    if (line[0] != '>')
    {
        LogError("ReadFasta: comment is not in correct format\n");
        exit(1);
    }

    line[strlen(line)-1] = '\0';
    if (comment != NULL)
        strcpy(comment, line + 1);

    int offset = 0;
    while (fgets(line, MaxLine, fp) != NULL)
    {
        if (line[0] == '>')
        {
            fseek(fp, -strlen(line), SEEK_CUR);
            break;
        }
        else
        {
            int len = strlen(line);
            while (isspace(line[len-1]))
                --len;
            copy(line, line + len, buf + offset);
            offset += len;
        }
    }

    seq->SetContent(buf, offset);
    return offset != 0;
}

void WriteFasta(FILE *fp, Sequence *seq, const char *comment)
{
    if (comment != NULL)
        fprintf(fp, ">%s\n", comment);

    int offset = 0;
    while (offset < seq->length)
    {
        int ch = seq->bases[offset + FastaMaxLine];
        seq->bases[offset + FastaMaxLine] = '\0';
        fprintf(fp, "%s\n", seq->bases + offset);
        seq->bases[offset + FastaMaxLine] = ch;
        offset += FastaMaxLine;
    }
}

uint64 CountSequences(std::FILE *fp, bool isCheck)
{
    uint64 index = 0;
    while (fgets(line, MaxLine, fp) != NULL)
    {
        if (line[0] == '>')
            ++index;
    }
    fseek(fp, 0, SEEK_SET);
    return index;
}

uint64 ReadHashGraph(std::FILE *fp, HashGraph *hashGraph, unsigned minCount)
{
    hashGraph->Clear();
    Sequence seq;

    uint64 index = 0;
    while (ReadFasta(fp, &seq, comment))
    {
        unsigned count = 0;
        for (int i = 0; comment[i]; ++i)
        {
            if (comment[i] == '_')
                comment[i] = ' ';
        }
        unsigned in, out;
        sscanf(comment, "%s %s %u %u %u", buf, buf, &count, &in, &out);
        seq.Encode();

        if (count >= minCount)
        {
            ++index;
            Kmer kmer = seq.GetKmer(0, seq.length);
            HashNode *node = hashGraph->InsertKmer(kmer);

            node->count = count;

            if (node->kmer != kmer)
                swap(in, out);

            node->in = in;
            node->out = out;
        }
    }

    return index;
}

uint64 WriteHashGraph(std::FILE *fp, HashGraph *hashGraph, unsigned minCount)
{
    uint64 index = 0;
    Sequence seq;

    for (unsigned i = 0; i < hashGraph->tableSize; ++i)
    {
        for (HashNode *node = hashGraph->table[i]; node != NULL; node = node->next)
        {
            if (!(node->status & HashGraphFlagDeadend) && node->count >= minCount)
            {
                Kmer kmer = node->kmer;
                seq.SetContent(kmer, kmerLength);
                seq.Decode();
                fprintf(fp, ">Kmer%lld_count_%u_%u_%u\n", index++, node->count, node->in, node->out);
                WriteFasta(fp, &seq);
            }
        }
    }

    return index;
}

// uint64 ReadContigGraph(std::FILE *fp, ContigGraph *contigGraph)
// {
//     contigGraph->numContigs = CountSequences(fp);
//     contigGraph->contigs = new Sequence[contigGraph->numContigs*2];
//     contigGraph->connections = new vector<Connection>[contigGraph->numContigs*2];
//     contigGraph->paths = new vector<Path>[contigGraph->numContigs*2];
//
//     Sequence seq;
//     for (unsigned i = 0; i < contigGraph->numContigs; ++i)
//     {
//         ReadFasta(fp, &seq);
//         seq.Encode();
//         contigGraph->contigs[i<<1] = seq;
//         seq.ReverseComplement();
//         contigGraph->contigs[(i<<1)|1] = seq;
//     }
//
//     return contigGraph->numContigs;
// }
