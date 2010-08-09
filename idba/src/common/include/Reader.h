#ifndef __READER_H_

#define __READER_H_

#include "globals.h"
#include "Read.h"

#include <string>
#include <fstream>
#include <cstdio>
#include <vector>

class Sequence;
class Contig;

class Reader
{
public:
    Reader(const std::string &filename, const std::string &mode);
    virtual ~Reader();

    virtual bool Read(Sequence &seq, std::string &comment) = 0;
    virtual int64 NumReads() = 0;

    bool Read(Contig &contig, std::string &comment);
    void Rewind() { std::fseek(fin, 0, SEEK_SET); }

protected:
    std::FILE *fin;
    char *line;
    char *buf;

private:
    Reader(const Reader &reader);
    const Reader &operator =(const Reader &reader);
};

class FastAReader: public Reader
{
public:
    FastAReader(const std::string &filename, const std::string &mode = "rb") 
        : Reader(filename, mode) {}
    ~FastAReader() {}

    virtual bool Read(Sequence &seq, std::string &comment);
    int64 NumReads();

private:
    FastAReader(const FastAReader &reader);
    const FastAReader &operator =(const FastAReader &reader);
};

class FastQReader: public Reader
{
public:
    FastQReader(const std::string &filename, const std::string &mode = "rb") 
        : Reader(filename, mode) {}
    ~FastQReader() {}

    virtual bool Read(Sequence &seq, std::string &comment);
    int64 NumReads();

private:
    FastQReader(const FastAReader &reader);
    const FastQReader &operator =(const FastAReader &reader);
};

#endif

