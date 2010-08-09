#include "globals.h"
#include "Sequence.h"
#include "Kmer.h"
#include "HashNode.h"
#include "HashTable.h"
#include "Log.h"
#include "Reader.h"
#include "Writer.h"

#include <omp.h>

#include <cstdio>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <queue>

using namespace std;

HashTable::Iterator::Iterator(HashTable *hashtable, int64 index, HashNode *current)
{ 
    this->hashtable = hashtable; 
    this->index = index;
    this->current = current;

    FindSolidPointer();
}

HashTable::Iterator &HashTable::Iterator::operator ++()
{
    if (hashtable == NULL)
        return *this;

    if (current != NULL)
        current = current->next;

    if (current == NULL)
    {
        ++index;
        FindSolidPointer();
    }

    return *this;
}

void HashTable::Iterator::FindSolidPointer()
{
    if (hashtable != NULL)
    {
        while (index < hashtable->table_size && hashtable->table[index] == NULL)
            ++index;
        
        if (index == hashtable->table_size)
            current = NULL;
        else
            current = hashtable->table[index];
    }
}

HashTable::HashTable(uint64 table_size, bool is_reversable)
{
    table = NULL;
    this->table_size = 0;
    this->is_reversable = is_reversable;
    num_nodes = 0;

    max_threads = omp_get_max_threads();
    if (max_threads > (int)MaxThreads)
    {
        max_threads = MaxThreads;
        omp_set_num_threads(MaxThreads);
    }

    mem_locks = new omp_lock_t[max_threads];
    heads = new HashNodePointer[max_threads];
    for (int i = 0; i < max_threads; ++i)
    {
        omp_init_lock(&mem_locks[i]);
        heads[i] = NULL;
    }

    omp_init_lock(&lock_alloc);
    for (uint64 i = 0; i < LockSize; ++i)
        omp_init_lock(&locks[i]);

    Reallocate(table_size);
}

HashTable::~HashTable()
{
    Clear();

    delete [] table;

    for (uint64 i = 0; i < LockSize; ++i)
        omp_destroy_lock(&locks[i]);

    omp_destroy_lock(&lock_alloc);
    for (int i = 0; i < max_threads; ++i)
        omp_destroy_lock(&mem_locks[i]);

    delete [] mem_locks;
    delete [] heads;
}


void HashTable::Clear()
{
    num_nodes = 0;

    for (unsigned i = 0; i < table_size; ++i)
        table[i] = NULL;
    
    reverse(backup.begin(), backup.end());
    for (unsigned i = 0; i < backup.size(); ++i)
    //for (unsigned i = backup.size()-1; i >= 0; --i)
        delete [] backup[i];

    for (int i = 0; i < max_threads; ++i)
        heads[i] = NULL;

    backup.resize(0);
}

void HashTable::ClearStatus()
{
#pragma omp parallel for
    for (int i = 0; i < table_size; ++i)
    {
        for (HashNode *node = table[i]; node; node = node->next)
            node->ClearStatus();
    }
}

void HashTable::SetDeadFlag()
{
#pragma omp parallel for
    for (int i = 0; i < table_size; ++i)
    {
        for (HashNode *node = table[i]; node; node = node->next)
            node->SetDeadFlag();
    }
}

void HashTable::ResetDeadFlag()
{
#pragma omp parallel for
    for (int i = 0; i < table_size; ++i)
    {
        for (HashNode *node = table[i]; node; node = node->next)
            node->ResetDeadFlag();
    }
}

void HashTable::SetUsedFlag()
{
#pragma omp parallel for
    for (int i = 0; i < table_size; ++i)
    {
        for (HashNode *node = table[i]; node; node = node->next)
            node->SetUsedFlag();
    }
}

void HashTable::ResetUsedFlag()
{
#pragma omp parallel for
    for (int i = 0; i < table_size; ++i)
    {
        for (HashNode *node = table[i]; node; node = node->next)
            node->ResetUsedFlag();
    }
}

KmerNode *HashTable::InsertKmer(const Kmer &kmer, int count)
{
    if (num_nodes > table_size)
    {
        Reallocate(table_size * 2);
    }

    Kmer key = kmer;

    if (is_reversable)
    {
        Kmer rev_comp = kmer;
        rev_comp.ReverseComplement();
        if (rev_comp < kmer)
        {
            key = rev_comp;
        }
    }

    uint64 index = key.Hash() % table_size;

    uint64 lockIndex = index & (LockSize - 1);
    omp_set_lock(&locks[lockIndex]);
    HashNode *p = table[index];

    while (p != NULL && p->kmer != key)
        p = p->next;

    if (p == NULL)
    {
        //p = NewNode(omp_get_thread_num());
        p = NewNode();
        p->SetContent(key);
        p->next = table[index];
        table[index] = p;

#pragma omp atomic
        num_nodes += 1;
    }

    p->Increase(count);
//    if (count + p->count < 65535)
//        p->count += count;
//    else
//        p->count = 65535;

    omp_unset_lock(&locks[lockIndex]);

    return p;
}

//HashNode *HashTable::GetByKmer(const Kmer &kmer)
//{
//    Kmer key = kmer;
//    
//    if (is_reversable)
//    {
//        Kmer rev_comp = kmer;
//        rev_comp.ReverseComplement();
//        if (rev_comp < kmer)
//            key = rev_comp;
//    }
//
//    uint64 index = key.Hash() % table_size;
//    HashNode *p = table[index];
//
//    while (p != NULL && p->kmer != key)
//        p = p->next;
//
//    return p;
//}

KmerNode *HashTable::GetNode(const Kmer &kmer)
{
    Kmer key = kmer;
    
    if (is_reversable)
    {
        Kmer rev_comp = kmer;
        rev_comp.ReverseComplement();
        if (rev_comp < kmer)
            key = rev_comp;
    }

    uint64 index = key.Hash() % table_size;
    HashNode *p = table[index];

    while (p != NULL && p->kmer != key)
        p = p->next;

    return p;
}

void HashTable::Save(std::vector<Kmer32> &kmers)
{
    kmers.resize(num_nodes);
    int index = 0;
    for (int i = 0; i < table_size; ++i)
    {
        for (HashNode *node = table[i]; node; node = node->next)
        {
            kmers[index++].data[0] = node->kmer.data[0];
        }
    }
}

void HashTable::Load(std::vector<Kmer32> &kmers)
{
//#pragma omp parallel for
    for (unsigned i = 0; i < kmers.size(); ++i)
    {
        Kmer kmer;
        kmer.data[0] = kmers[i].data[0];
        InsertKmer(kmer, 0);
    }
}

void HashTable::Save(Writer &writer)
{
//    for (int i = 0; i < table_size; ++i)
//    {
//        for (HashNode *node = table[i]; node; node = node->next)
//        {
//            Sequence seq = node->kmer;
//            seq.Decode();
//            writer.WriteFormat(seq, "kmer_%d_%u", node->Count(), node->edges);
//        }
//    }
}

void HashTable::Load(Reader &reader)
{
//    Sequence seq;
//    string comment;
//    while (reader.Read(seq, comment))
//    {
//        seq.Encode();
//        Kmer kmer = seq.GetKmer(0, kmerLength);
//        HashNode *node = InsertKmer(kmer);
//        int count;
//        unsigned edges;
//        for (unsigned i = 0; i < comment.size(); ++i)
//        {
//            if (comment[i] == '_')
//                comment[i] = ' ';
//        }
//        char tmp[100];
//        sscanf(comment.c_str(), "%s %d %u", tmp, &count, &edges);
//        node->count = count;
//        node->edges = edges;
////        node->in = in;
////        node->out = out;
////        node->count = count;
//    }
}

HashNode *HashTable::NewNode(int index)
{
    omp_set_lock(&mem_locks[index]);
    if (heads[index] == NULL)
    {
        heads[index] = AllocateNodes();
    }
    HashNode *node = heads[index];
    heads[index] = node->next;
    omp_unset_lock(&mem_locks[index]);
    return node;
}

void HashTable::FreeNode(HashNode *node, int index)
{
    omp_set_lock(&mem_locks[index]);
    node->next = heads[index];
    heads[index] = node;
    omp_unset_lock(&mem_locks[index]);
}

HashNode *HashTable::AllocateNodes()
{
    omp_set_lock(&lock_alloc);
    HashNode *buffer = new HashNode[HashNodeBufferSize];
    for (unsigned i = 0; i+1 < HashNodeBufferSize; ++i)
    {
        buffer[i].next = &buffer[i+1];
    }
    buffer[HashNodeBufferSize-1].next = NULL;
    backup.push_back(buffer);
    omp_unset_lock(&lock_alloc);
    return buffer;
}

void HashTable::Reallocate(int64 new_table_size)
{
    omp_set_lock(&locks[0]);

    if (table_size == new_table_size)
    {
        omp_unset_lock(&locks[0]);
        return;
    }

    for (unsigned i = 1; i < LockSize; ++i)
        omp_set_lock(&locks[i]);

    HashNode **new_table = new HashNodePointer[new_table_size];

    if (new_table == NULL)
    {
        LogError("HashGraph::Reallocate() not enough memory\n");
        exit(1);
    }

    for (int64 i = 0; i < new_table_size; ++i)
        new_table[i] = NULL;

    for (int64 i = 0; i < table_size; ++i)
    {
        HashNode *node = table[i];
        while (node != NULL)
        {
            HashNode *next = node->next;
            uint64 index = node->kmer.Hash() % new_table_size;
            node->next = new_table[index];
            new_table[index] = node;
            node = next;
        }
    }

    delete [] table;

    table = new_table;
    table_size = new_table_size;

    for (unsigned i = 0; i < LockSize; ++i)
        omp_unset_lock(&locks[i]);
}
