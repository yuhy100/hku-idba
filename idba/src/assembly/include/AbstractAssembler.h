#ifndef __ABSTRACT_ASSEMBLER_H_

#define __ABSTRACT_ASSEMBLER_H_

#include "globals.h"
#include "Option.h"
#include "Read.h"
#include "Reader.h"
#include "Sequence.h"

#include <vector>

class Option;
class Contig;
class Read;

class AbstractAssembler
{
public:
    AbstractAssembler(const Option &option, int64 read_buffer = 0);
    virtual ~AbstractAssembler() {}

    std::vector<Sequence> &GetLongReads()
    { return long_reads; }

    std::vector<Read> &GetShortReads();
    void RewindShortReads()
    { read_reader->Rewind(); read_num = 0; }

    virtual void Assemble(std::vector<Contig> &contigs) = 0;

protected:
    Option option;

private:
    void ReadBunchOfReads(std::vector<Read> &reads, int64 expected_num);

    std::vector<Sequence> long_reads;
    std::vector<Read> reads;
    std::vector<Read> empty_reads;
    int read_buffer;
    Reader *read_reader;
    int64 read_num;
};

#endif
