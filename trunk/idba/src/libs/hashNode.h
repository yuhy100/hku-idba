#ifndef __HASHNODE_H_

#define __HASHNODE_H_

#include "globals.h"
#include "kmer.h"

#include <algorithm>

struct HashNode
{
    HashNode *next;
    Kmer kmer;
    unsigned short status;
    unsigned short count;
    int data[0];
    unsigned char inDepth;
    unsigned char outDepth;
    unsigned char in;
    unsigned char out;

    bool GetStatus(unsigned flag)
    { return status & flag; }
    void SetStatus(unsigned flag)
    { status |= flag; }
    void ResetStatus(unsigned flag)
    { status &= ~flag; }
    void ToggleStatus(unsigned flag)
    { status ^= flag; }

    int &Data() { return data[0]; }

    void SetContent(const Kmer &kmer)
    {
        next = NULL;
        this->kmer = kmer;
        status = 0;
        inDepth = outDepth = 0;
        count = 0;
        in = out = 0;
    }

    static HashNode *start;
    static unsigned totalNodes;

    static HashNode *NewNode()
    {
	    if (start == NULL)
	    {
            AddBuffer();
	    }

	    HashNode *p = start;
	    start = start->next;
	    return p;
    }

    static void FreeNode(HashNode *node)
    {
	    node->next = start;
	    start = node;
    }

    static void AddBuffer();
};

struct HashNodeAdapter
{
    HashNode *node;
    bool isReverse;

    HashNodeAdapter(HashNode *node) { SetHashNode(node); }
    HashNodeAdapter(HashNode *node, const Kmer &kmer) { SetHashNode(node, kmer); }

    void SetHashNode(HashNode *node)
    {
        this->node = node;
        isReverse = false;
    }

    void SetHashNode(HashNode *node, const Kmer &kmer)
    {
        SetHashNode(node);
        isReverse = (kmer != node->kmer);
    }

    void SetDirection(const Kmer &kmer) { isReverse = (kmer != node->kmer); }
    void ReverseComplement() { isReverse = !isReverse; }

    unsigned char &OutEdges() { return (!isReverse ? node->out : node->in); }
    unsigned char &InEdges() { return (!isReverse ? node->in : node->out); }
    unsigned char &OutDepth() { return (!isReverse ? node->outDepth : node->inDepth); }
    unsigned char &InDepth() { return (!isReverse ? node->inDepth : node->outDepth); }

    void AddOutEdge(unsigned x) { OutEdges() |= (1 << x); }
    void RemoveOutEdge(unsigned x) { OutEdges() &= ~(1 << x); }
    void AddInEdge(unsigned x) { InEdges() |= (1 << x); }
    void RemoveInEdge(unsigned x) { InEdges() &= ~(1 << x); }

    int OutDegree();
    int InDegree();
};

#endif
