#ifndef __HASH_ALIGN_H_

#define __HASH_ALIGN_H_

#include "HashTable.h"
#include "Contig.h"
#include "ContigNode.h"

struct Alignment
{
    int readLength;
    int readOffset;
    int contigId;
    int contigOffset;
    int contigLength;
    int length;
    bool isReverse;

    void ReverseComplement()
    {
        readOffset = readLength - (readOffset + length);
        contigOffset = contigLength - (contigOffset + length);
        //contigId ^= 1;
        isReverse = !isReverse;
    }
};

struct Position
{
    int id;
    int offset;
    int length;
    bool isReverse;

    void ReverseComplement()
    {
        offset = length - (offset + kmerLength);
        isReverse = !isReverse;
    }
};

class HashAlign: public HashTable
{
public:
    HashAlign(int64 tableSize = MaxHashTable): HashTable(tableSize, true) {}

    void Initialize(std::vector<Contig> &contigs);
    void Initialize(std::vector<ContigNode> &contigs);
    void InsertSequence(const Sequence &seq);
    void InsertSequence(const Sequence *seq, int id);
    bool AlignKmer(const Kmer &kmer, Position &pos);
    void AlignSequence(const Sequence *seq, std::vector<Alignment> &alignments);

private:
    HashAlign(const HashAlign &);
    const HashAlign &operator =(const HashAlign &);

    void BuildTable(std::vector<Contig> &contigs);
    void BuildPositions(std::vector<Contig> &contigs);

    std::vector<Position> positions;
};

#endif
