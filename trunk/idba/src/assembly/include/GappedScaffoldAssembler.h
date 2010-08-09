#ifndef __GAPPED_SCAFFOLD_ASSEMBLER_H_

#define __GAPPED_SCAFFOLD_ASSEMBLER_H_

#include "globals.h"
#include "AbstractAssembler.h"
#include "ScaffoldGraph.h"
#include "HashAlign.h"

class Contig;

class GappedScaffoldAssembler: public AbstractAssembler
{
public:
    GappedScaffoldAssembler(const Option &option, int read_buffer = 0)
        : AbstractAssembler(option, read_buffer) {}

    virtual void Assemble(std::vector<Contig> &contigs);

private:
    ScaffoldGraph scaffold_graph;
    HashAlign aligner;
};

#endif
