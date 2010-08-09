#ifndef __SIMPLE_ASSMBLER_H_

#define __SIMPLE_ASSMBLER_H_

#include "globals.h"
#include "HashGraph.h"
#include "AbstractAssembler.h"

#include <vector>

class Contig;

class SimpleAssembler: public AbstractAssembler
{
public:
    SimpleAssembler(const Option &option, int64 read_buffer = 0)
        : AbstractAssembler(option, read_buffer) 
    { hash_graph = new HashGraph; }

    virtual ~SimpleAssembler() 
    { delete hash_graph; }

    virtual void Assemble(std::vector<Contig> &contigs);
    
private:
    void BuildKmerNodes();
    void AddInternalKmerNodes();

    HashGraph *hash_graph;
};

#endif
