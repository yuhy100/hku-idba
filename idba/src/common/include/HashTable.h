#ifndef __HASH_TABLE_H_

#define __HASH_TABLE_H_

#include "globals.h"
#include "KmerNode.h"
#include "HashNode.h"

#include <omp.h>

#include <vector>

class Sequence;
class Kmer;
class Reader;
class Writer;

#include <iostream>

class HashTable
{
    friend class SimpleAssembler;

public:
    typedef HashNode *HashNodePointer;

    class Iterator
    {
    public:
        Iterator(HashTable *hashtable = NULL, int64 index = 0, HashNode *current = NULL);

        Iterator &operator ++();

        KmerNode *operator *() { return current; }
        bool operator == (const Iterator &iter)
        {
            return hashtable == iter.hashtable && index == iter.index && current == iter.current;
        }
        bool operator != (const Iterator &iter)
        {
            return !((*this) == iter);
        }

    private:
        void FindSolidPointer();

        HashTable *hashtable;
        int64 index;
        HashNode *current;
    };

    HashTable(uint64 table_size = MaxHashTable, bool is_reversable = false);
    ~HashTable();

    Iterator Begin() { return Iterator(this, 0, 0); }
    Iterator End() { return Iterator(this, table_size, 0); }

    void Clear();
    void ClearStatus();
    void SetDeadFlag();
    void ResetDeadFlag();
    void SetUsedFlag();
    void ResetUsedFlag();

    KmerNode *InsertKmer(const Kmer &kmer, int count = 1);
    KmerNode *GetNode(const Kmer &kmer);

    int64 NumNodes() const { return num_nodes; }

    void Save(std::vector<Kmer32> &kmers);
    void Load(std::vector<Kmer32> &kmers);

    void Save(Writer &writer);
    void Load(Reader &reader);

    static const uint64 HashNodeBufferSize = (1 << 12);
    static const uint64 MaxHashTable = 999997;
    static const uint64 LockSize = (1 << 20);
    static const uint64 MaxThreads = (1 << 6);

protected:
    HashNode *NewNode() { return NewNode(omp_get_thread_num()); }
    void FreeNode(HashNode *node) { FreeNode(node, omp_get_thread_num()); }

    HashNode *NewNode(int index);
    void FreeNode(HashNode *node, int index);

    HashNode **table;
    int64 table_size;
    int64 num_nodes;

private:
    HashTable(const HashTable &);
    const HashTable &operator =(const HashTable &);

    HashNode *AllocateNodes();
    void Reallocate(int64 new_table_size);

    bool is_reversable;
    int max_threads;
    omp_lock_t locks[LockSize];
    HashNode **heads;
    omp_lock_t *mem_locks;
    std::vector<HashNode *> backup;
    omp_lock_t lock_alloc;

};

#endif

