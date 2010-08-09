#ifndef __ITERATIVE_ASSEMBLER_H_

#define __ITERATIVE_ASSEMBLER_H_

#include "globals.h"
#include "AbstractAssembler.h"
#include "ContigGraph.h"

#include <vector>

class Contig;

class IterativeAssembler: public AbstractAssembler
{
public:
    IterativeAssembler(const Option &option, int64 read_buffer = 0)
        : AbstractAssembler(option, read_buffer) {}

    virtual void Assemble(std::vector<Contig> &contigs);

private:
    void Iterate();

    ContigGraph contig_graph;
    std::vector<Contig> tmp_contigs;
};

#endif
