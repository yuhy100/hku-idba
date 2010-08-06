#include "globals.h"
#include "hashNode.h"

using namespace std;

HashNode *HashNode::start = NULL;
unsigned HashNode::totalNodes = 0;

static unsigned bits[] =
{
    0, 1, 1, 2, 1, 2, 2, 3,
    1, 2, 2, 3, 2, 3, 3, 4,
};

void HashNode::AddBuffer()
{
    totalNodes += HashNodeBufferSize;
	HashNode *nodes = new HashNode[HashNodeBufferSize];
	for (unsigned i = 0; i < HashNodeBufferSize-1; ++i)
		nodes[i].next = &nodes[i+1];
	nodes[HashNodeBufferSize-1].next = start;
	start = &nodes[0];
}

int HashNodeAdapter::OutDegree()
{
    return bits[OutEdges()];
}

int HashNodeAdapter::InDegree()
{
    return bits[InEdges()];
}