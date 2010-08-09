#ifndef __WRITER_H_

#define __WRITER_H_

#include "globals.h"

#include <string>
#include <fstream>

class Sequence;
class Contig;

class Writer
{
public:
    Writer(const std::string &filename, const std::string &mode);
    virtual ~Writer();

    virtual bool Write(const Sequence &seq, const std::string &comment) = 0;
    bool WriteFormat(const Sequence &seq, const char *fmt, ...);
    bool Write(const Contig &contig, const std::string &comment);
    bool WriteFormat(const Contig &contig, const char *fmt, ...);

protected:
    std::FILE *fout;
    char *line;

private:
    Writer(const Writer &writer);
    const Writer &operator =(const Writer &writer);
};

class FastAWriter: public Writer
{
public:
    FastAWriter(const std::string &filename, const std::string &mode = "wb")
        : Writer(filename, mode) {}
    ~FastAWriter() {}

    virtual bool Write(const Sequence &seq, const std::string &comment);

private:
    FastAWriter(const FastAWriter &writer);
    const FastAWriter &operator =(const FastAWriter &writer);
};

class FastQWriter: public Writer
{
public:
    FastQWriter(const std::string &filename, const std::string &mode = "wb")
        : Writer(filename, mode) {}
    ~FastQWriter() {}

    virtual bool Write(const Sequence &seq, const std::string &comment);
    bool Write(const Sequence &seq, const std::string &comment, int id);

private:
    FastQWriter(const FastAWriter &writer);
    const FastQWriter &operator =(const FastAWriter &writer);
};

#endif

