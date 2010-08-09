#include "globals.h"
#include "Sequence.h"
#include "Contig.h"
#include "HashNode.h"
#include "Reader.h"
#include "Utils.h"
#include "Writer.h"
#include "HashGraph.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <fstream>


using namespace std;

char line[10000000];

int main(int argc, char *argv[])
{
    cout << sizeof(AbstractNode) << endl;
    cout << sizeof(KmerNode) << endl;
    cout << sizeof(HashNode) << endl;

    HashNode node;
    node.Data() = 10;
    cout << node.Data() << endl;

    node.Data() = 65537;
    cout << node.Data() << endl;

    node.Data() = 65536;
    cout <<  node.Data() << endl;

    node.SetInEdges(3);
    node.SetOutEdges(9);

    node.RemoveInEdge(1);

    HashGraph *hash_graph = new HashGraph;
    FastAReader reader(argv[1]);
    Sequence seq;
    string comment;
    reader.Read(seq, comment);
    seq.Encode();

    hash_graph->InsertSequence(seq);
    cout << hash_graph->NumNodes() << endl;

//    vector<Kmer32> kmers;
//    hash_graph->Save(kmers);


    hash_graph->Clear();
    cout << hash_graph->NumNodes() << endl;
    delete hash_graph;


    while (true)
        ;



//    cout << node.InEdges() << " " << node.OutEdges() << endl;
//    node.ReverseComplement();
//    cout << node.InEdges() << " " << node.OutEdges() << endl;
//    node.ReverseComplement();
//    cout << node.InEdges() << " " << node.OutEdges() << endl;


//    kmerLength = 50;
//
//    vector<Contig> contigs;
//    FastAReader reader("tmp.fa");
//    HashGraph hashGraph;
//    hashGraph.Load(reader);
//    hashGraph.Assemble(contigs);
//
//    FastAWriter writer("out.fa");
//    for (int i = 0; i < contigs.size(); ++i)
//    {
//        contigs[i].Decode();
//        writer.WriteFormat(contigs[i], "contig%d_%d", i, contigs[i].Size());
//    }

//    FastAReader reader(argv[1]);
//    FastAWriter writer(argv[2]);
//    
//    //FILE *fp = OpenFile(argv[1], "rb");
//    FILE *fp = OpenFile(argv[2], "wb");
//    ofstream fout(argv[2], ios_base::binary | ios_base::out);
//
//    Sequence seq;
//    string comment;
//    int count = 0;
//    while (reader.Read(seq, comment))
//    {
//        //writer.Write(seq, "%s_%d", comment.c_str(), count);
//        writer.WriteFormat(seq, "%s_%d", comment.c_str(), count);
//        //Writer *wr = &writer;
//        //wr->Write(seq, comment, "_%d", count);
//        //writer.Write(seq, comment, count);
//        //WriteFastA(fp, &seq, comment.c_str());
//        //fout << seq << endl;
//        //fprintf(fp, "%s", seq.ToCString());
//        ++count;
//    }
//    cout << count << endl;
//    //cout << reader.Size() << endl;
//
//    return 0;
}

