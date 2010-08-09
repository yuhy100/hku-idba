#include "globals.h"
#include "Sequence.h"
#include "Contig.h"
#include "Reader.h"

#include <cstdio>
#include <string>
#include <iostream>
#include <algorithm>
#include <cctype>

using namespace std;

Reader::Reader(const string &filename, const string &mode)
{
    fin = fopen(filename.c_str(), mode.c_str());
    if (fin == NULL)
    {
        cerr << "can't open " << filename << endl;
        throw exception();
    }
    line = new char[MaxLine];
    buf = new char[MaxLine];
}

Reader::~Reader()
{
    if (fin != NULL)
        fclose(fin);
    delete [] line;
    delete [] buf;
}

bool Reader::Read(Contig &contig, string &comment)
{
    if (!Read((Sequence &)contig, comment))
        return false;

    int index = comment.size()-1;
    while (isdigit(comment[index]))
        --index;
    contig.sum_coverage = atoi(comment.c_str() + index + 1);
    return true;
}

bool FastAReader::Read(Sequence &seq, string &comment)
{
    if (fgets(line, MaxLine, fin) == NULL)
        return false;

    if (line[0] != '>')
    {
        cerr << "ReadFastA: comment is not in correct format" << endl;
        throw exception();
    }

    line[strlen(line)-1] = '\0';
    comment = line + 1;
    seq.Clear();

    int offset = 0;
    while (fgets(line, MaxLine, fin) != NULL)
    {
        if (line[0] == '>')
        {
            fseek(fin, -strlen(line), SEEK_CUR);
            break;
        }

        int len = strlen(line);
        while (isspace(line[len-1]))
            --len;
        line[len] = '\0';

        if (len + offset + 1 > MaxLine)
        {
            seq += buf;
            offset = 0;
        }

        strcpy(buf + offset, line);
        offset += len;
    }

    if (offset > 0)
    {
        seq += buf;
        offset = 0;
    }

    return true;
}

int64 FastAReader::NumReads()
{
    int64 index = ftell(fin);
    fseek(fin, 0, SEEK_SET);

    uint64 num = 0;
    while (fgets(line, MaxLine, fin) != NULL)
    {
        if (line[0] == '>')
            ++num;
    }

    fseek(fin, index, SEEK_SET);
    return num;
}

bool FastQReader::Read(Sequence &seq, string &comment)
{
    if (fgets(line, MaxLine, fin) == NULL)
        return false;

    if (line[0] != '@')
    {
        cerr << "ReadFastQ: comment is not in correct format" << endl;
        throw exception();
    }

    line[strlen(line)-1] = '\0';
    comment = line + 1;
    seq.Clear();

    int offset = 0;
    while (fgets(line, MaxLine, fin) != NULL)
    {
        if (line[0] == '+')
        {
            fgets(line, MaxLine, fin);
            line[strlen(line)-1] = '\0';
            break;
        }

        int len = strlen(line);
        while (isspace(line[len-1]))
            --len;
        line[len] = '\0';

        if (len + offset + 1 > MaxLine)
        {
            seq += buf;
            offset = 0;
        }

        strcpy(buf + offset, line);
        offset += len;
    }

    if (offset > 0)
    {
        seq += buf;
        offset = 0;
    }

    return true;
}

int64 FastQReader::NumReads()
{
    int64 index = ftell(fin);
    fseek(fin, 0, SEEK_SET);

    uint64 num = 0;
    while (fgets(line, MaxLine, fin) != NULL)
    {
        if (line[0] == '@')
            ++num;
    }

    fseek(fin, index, SEEK_SET);
    return num;
}

