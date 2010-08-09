#include "globals.h"
#include "Sequence.h"
#include "Contig.h"
#include "Writer.h"

#include <cstdio>
#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cstdarg>

using namespace std;

Writer::Writer(const string &filename, const string &mode)
{
    fout = fopen(filename.c_str(), mode.c_str());
    if (fout == NULL)
    {
        cerr << "can't open " << filename << endl;
        throw exception();
    }
    line = new char[MaxLine];
}

Writer::~Writer()
{
    if (fout != NULL)
        fclose(fout);
    delete [] line;
}

bool Writer::WriteFormat(const Sequence &seq, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsprintf(line, fmt, ap);
    va_end(ap);

    return Write(seq, line);
}

bool Writer::Write(const Contig &contig, const string &comment)
{
    return WriteFormat(contig, "%s", comment.c_str());
}

bool Writer::WriteFormat(const Contig &contig, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsprintf(line, fmt, ap);
    va_end(ap);

    sprintf(line + strlen(line), "_%d", contig.sum_coverage);

    return Write((const Sequence &)contig, line);
}


bool FastAWriter::Write(const Sequence &seq, const string &comment)
{
    fprintf(fout, ">%s\n", comment.c_str());
    fprintf(fout, "%s\n", seq.ToCString());

    return true;
}

bool FastQWriter::Write(const Sequence &seq, const string &comment)
{
    fprintf(fout, "@%s\n", comment.c_str());
    fprintf(fout, "%s\n", seq.ToCString());
    fprintf(fout, "+\n");
    string quality;
    quality.resize(seq.Size());
    fill(quality.begin(), quality.end(), 'a');
    fprintf(fout, "%s\n", quality.c_str());

    return true;
}

