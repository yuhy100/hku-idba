#ifndef __HASH_NODE_H_

#define __HASH_NODE_H_

#include "globals.h"
#include "Kmer.h"
#include "BitOperation.h"
#include "AbstractNode.h"
#include "KmerNode.h"

#include <algorithm>


class HashNode: public KmerNode
{
public:
    friend class HashTable;
    friend class HashGraph;
    HashNode () { next = NULL; }

    void SetContent(const Kmer &kmer)
    {
        KmerNode::SetContent(kmer);
        next = NULL;
    }

    HashNode *GetNext() { return next; }

private:

    HashNode *next;
};

#endif
