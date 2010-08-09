#include "globals.h"
#include "Kmer.h"
#include "Sequence.h"
#include "Utils.h"
#include "HashNode.h"
#include "HashGraph.h"
#include "ContigGraph.h"
#include "ContigNode.h"
#include "Log.h"
#include "ContigBranchGroup.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <vector>
#include <cstring>

using namespace std;


void ContigGraph::Refresh()
{
    hash_graph->SetDeadFlag();
#pragma omp parallel for
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        //contigs[i].ResetUsedFlag();
        nodes[i].ClearLock();

        KmerNode *p = hash_graph->GetNode(nodes[i].GetBeginKmer());
        if (p == NULL)
        {
            cout << "Error" << endl;
            Contig contig = nodes[i].GetContig();
            contig.Decode();
            cout << contig << endl;
            contig.ReverseComplement();
            cout << contig << endl;
        }
        p->Data() = i;
        p->ResetDeadFlag();

        p = hash_graph->GetNode(nodes[i].GetEndKmer());
        if (p == NULL)
        {
            cout << "Error" << endl;
            Contig contig = nodes[i].GetContig();
            contig.Decode();
            cout << contig << endl;
            contig.ReverseComplement();
            cout << contig << endl;
        }
        p->Data() = i;
        p->ResetDeadFlag();
    }
    hash_graph->Refresh();

#pragma omp parallel for
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        ContigNodeAdapter contigAdp(&nodes[i]);
        for (int strand = 0; strand < 2; ++strand)
        {
            KmerNodeAdapter adp = hash_graph->GetNodeAdapter(contigAdp.GetEndKmer());
            contigAdp.SetOutEdges(adp.OutEdges());
           

            for (int x = 0; x < 4; ++x)
            {
                if ((1 << x) & adp.OutEdges())
                {
                    ContigNodeAdapter neighbor = GetNeighbor(contigAdp, x);
                    if (neighbor.IsNull())
                    {
                        adp.RemoveOutEdge(x);
                        contigAdp.RemoveOutEdge(x);

#ifdef DEBUG
                        //LogMessage("occational hit\n");
#endif
                    }
                }
            }
            
            contigAdp.ReverseComplement();
        }
    }
}

//bool ContigGraph::Check()
//{
//#pragma omp parallel for
//    for (int64 i = 0; i < (int64)contigs.size(); ++i)
//    {
//        ContigNodeAdapter contigAdp(&contigs[i]);
//        if (contigs[i].IsDead())
//            continue;
//
//        for (int strand = 0; strand < 2; ++strand)
//        {
//            HashNode *p = hash_graph->GetByKmer(contigAdp.outKmer);
//            KmerNodeAdapter adp(p, contigAdp.outKmer);
//            contigAdp.out = adp.OutEdges();
//
//            if (contigAdp.outKmer != adp.GetKmer())
//            {
//                cout << "vertex doesn't match the contig" << endl;
//                cout << contigAdp.outKmer << endl;
//                cout << adp.GetKmer() << endl;
//                exit(1);
//            }
//
//            if (p->Data() != i)
//            {
//                int j = p->Data();
//                cout << i << " " << j << endl;
//                contigs[i].Decode();
//                contigs[j].Decode();
//                cout << contigs[i] << endl;
//                cout << contigs[j] << endl;
//
//                exit(1);
//            }
//
//            for (int j = 0; j < 4; ++j)
//            {
//                if ((1 << j) & contigAdp.out)
//                {
//                    Kmer x = contigAdp.outKmer;
//                    x.AddRight(j);
//
//                    HashNode *q = hash_graph->GetByKmer(x);
//                    if (q == NULL)
//                        continue;
//
//                    ContigNodeAdapter to(&contigs[hash_graph->GetByKmer(x)->Data()]);
//                    if (to.inKmer != x)
//                    {
//                        to.ReverseComplement();
//                        if (to.inKmer != x)
//                        {
//                            cout << "kmer " << kmerLength << endl;
//                            cout << "occational hit" << endl;
//                            cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
//                            cout << x << " " << to.inKmer << endl;
//                            to.ReverseComplement();
//                            cout << x << " " << to.inKmer << endl;
//                            to.ReverseComplement();
//
//                            ContigNode c = contigs[i];
//                            if (c.outKmer != contigAdp.outKmer)
//                            {
//                                c.ReverseComplement();
//                            }
//                            c.Decode();
//                            cout << c << endl;
//                            c = *to.contig;
//                            c.Decode();
//                            cout << c << endl;
//                            c.ReverseComplement();
//                            cout << c << endl;
//                        }
//                    }
//                }
//            }
//
//            contigAdp.ReverseComplement();
//        }
//    }
//
//    return true;
//}

double ContigGraph::AverageCoverage()
{
    double sum = 0;
    int64 count = 0;
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        sum += nodes[i].GetContig().sum_coverage;
        count += nodes[i].GetContig().Size() - kmerLength + 1;
    }
    return sum / count;
}

int ContigGraph::RemoveDeadEnd(int minLength)
{
    int deadend = 0;
#pragma omp parallel for
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        if (nodes[i].GetContig().Size() < kmerLength + minLength - 1 && 
                (nodes[i].InDegree() == 0 || nodes[i].OutDegree() == 0))
        {
            nodes[i].SetDeadFlag();

#pragma omp atomic
            ++deadend;
        }
    }

    Refresh();

    return deadend;
}

int ContigGraph::RemoveLowCoverageContigs(double c)
{
    int remove = 0;
#pragma omp parallel for
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        if (nodes[i].GetContig().Coverage() < c)
        {
            nodes[i].SetDeadFlag();

#pragma omp atomic
            ++remove;
        }
    }

    Refresh();

    return remove;
}

int ContigGraph::RemoveBubble()
{
    int bubbles = 0;
    for (unsigned k = 0; k < nodes.size(); ++k)
    {
        if (nodes[k].IsDead())
            continue;

        ContigNodeAdapter begin(&nodes[k]);
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigBranchGroup branch_group(this, begin);
            if (branch_group.Search())
            {
                branch_group.Merge();
                ++bubbles;
            }

            begin.ReverseComplement();
        }
    }
    
    Refresh();

    return bubbles;
}

void ContigGraph::MergeContigs()
{
    vector<Contig> new_contigs;
    Assemble(new_contigs);
    SetContigs(new_contigs);
    Refresh();
}

int ContigGraph::Assemble(vector<Contig> &result_contigs)
{
    result_contigs.resize(0);

    omp_lock_t lockContigs;
    omp_init_lock(&lockContigs);

#pragma omp parallel for
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        if (nodes[i].GetContig().Size() == kmerLength 
                && nodes[i].GetBeginKmer().IsPalindrome()
                && !nodes[i].IsDead())
        {
            nodes[i].SetDeadFlag();

            omp_set_lock(&lockContigs);
            result_contigs.push_back(nodes[i].GetContig());
            omp_unset_lock(&lockContigs);
        }
    }

    int total = 0;
#pragma omp parallel for
    for (int i = 0; i < (int64)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        if (!nodes[i].Lock(omp_get_thread_num()))
            continue;

        ++total;
        ContigNodeAdapter current(&nodes[i]);
        ContigPath path;
        path.Append(current);

        Contig contig;
        for (int strand = 0; strand < 2; ++strand)
        {
            while (true)
            {
                current = path.GetEndNodeAdapter();
                ContigNodeAdapter next;

                if (!GetNextNodeAdapter(current, next))
                    break;

                if (next.IsDead())
                    break;

                if (next.GetLockID() == omp_get_thread_num() && IsLoop(path, next))
                    break;

                if (!next.LockPreempt(omp_get_thread_num()))
                    goto FAIL;

                path.Append(next);
            }

            path.ReverseComplement();
        }

        path.Merge(contig);

        omp_set_lock(&lockContigs);
        result_contigs.push_back(contig);
        omp_unset_lock(&lockContigs);

FAIL:
        ;
    }

    omp_destroy_lock(&lockContigs);

    return result_contigs.size();
}

int ContigGraph::Assemble(vector<Contig> &result_contigs, vector<Kmer> &branches)
{
    result_contigs.resize(0);
    branches.resize(0);

    omp_lock_t lockKmers;
    omp_init_lock(&lockKmers);

    int64 size = nodes.size();
#pragma omp parallel for
    for (int64 i = 0; i < size; ++i)
    {
        if (!nodes[i].IsDead() && nodes[i].GetSize() >= kmerLength)
        {
            ContigNodeAdapter contig_adp(&nodes[i]);
            for (int strand = 0; strand < 2; ++strand)
            {
                Kmer kmer = contig_adp.GetEndKmer();
                unsigned edges = contig_adp.OutEdges();
                
                for (int j = 0; j < 4; ++j)
                {
                    if (edges & (1 << j))
                    {
                        Kmer x = kmer;
                        if (x.IsPalindrome() && contig_adp.GetSize() > kmerLength 
                                && contig_adp.GetNucleotide(contig_adp.GetSize() - kmerLength - 1) == (3U - j))
                            continue;

                        x.SetBase(kmerLength, j);
                        omp_set_lock(&lockKmers);
                        branches.push_back(x);
                        omp_unset_lock(&lockKmers);
                    }
                }

                contig_adp.ReverseComplement();
            }
        }
    }

    ++kmerLength;

#pragma omp parallel for
    for (int64 i = 0; i < (int64)branches.size(); ++i)
    {
        Kmer revComp = branches[i];
        revComp.ReverseComplement();
        if (revComp < branches[i])
            branches[i] = revComp;
    }

    sort(branches.begin(), branches.end());
    branches.erase(unique(branches.begin(), branches.end()), branches.end());

    --kmerLength;

    return Assemble(result_contigs);
}

ContigNodeAdapter ContigGraph::GetNeighbor(const ContigNodeAdapter &node, int x)
{
    Kmer kmer = node.GetEndKmer();
    kmer.AddRight(x);

    KmerNode *p = hash_graph->GetNode(kmer);
    if (p == NULL)
        return ContigNodeAdapter(NULL);

    ContigNodeAdapter neighbor(&nodes[p->Data()]);
    if (neighbor.GetBeginKmer() != kmer)
        neighbor.ReverseComplement();

    if (neighbor.GetBeginKmer() != kmer)
        return ContigNodeAdapter(NULL);
    else
        return neighbor;
} 

void ContigGraph::GetNeighbors(const ContigNodeAdapter &node, 
        vector<ContigNodeAdapter> &neighbors)
{
    neighbors.resize(0);
    for (int x = 0; x < 4; ++x)
    {
        if ((1 << x) & node.OutEdges())
        {
            Kmer kmer = node.GetEndKmer();
            kmer.AddRight(x);

            int j = hash_graph->GetNode(kmer)->Data();
            ContigNodeAdapter neighbor(&nodes[j]);
            if (neighbor.GetBeginKmer() != kmer)
                neighbor.ReverseComplement();

            neighbors.push_back(neighbor);
        }
    }
}

//void ContigGraph::RemoveContigNode(ContigNode &node)
//{
//    HashNode *p = hash_graph->GetNode(node.GetBeginKmer());
//    p->SetDeadFlag();
//
//    p = hash_graph->GetNode(node.GetEndKmer());
//    p->SetDeadFlag();
//
//    node.SetDeadFlag();
//}
//

void ContigGraph::BuildVertices()
{
    hash_graph->Clear();

#pragma omp parallel for
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        if (nodes[i].GetContig().Size() < kmerLength)
        {
            nodes[i].SetDeadFlag();
            continue;
        }

        nodes[i].Clear();
        hash_graph->InsertKmer(nodes[i].GetBeginKmer());
        hash_graph->InsertKmer(nodes[i].GetEndKmer());
    }
}

void ContigGraph::SetContigs(vector<Contig> &contigs)
{
    this->nodes.resize(contigs.size());
#pragma omp parallel for
    for (int i = 0; i < (int)contigs.size(); ++i)
    {
        this->nodes[i].SetContent(contigs[i]);
        this->nodes[i].Data() = i;
    }
}

void ContigGraph::AddBranches(vector<Kmer> &branches)
{
    int old_size = nodes.size();
    nodes.resize(nodes.size() + branches.size());
#pragma omp parallel for
    for (int i = 0; i < (int)branches.size(); ++i)
    {
        Contig contig(branches[i]);
        nodes[old_size + i].SetContent(contig);
        nodes[old_size + i].Data() = old_size + i;
    }
}

bool ContigGraph::GetNextNodeAdapter(ContigNodeAdapter &current, ContigNodeAdapter &next)
{
    if (current.OutDegree() != 1)
        return false;

    Kmer kmer = current.GetEndKmer();
    kmer.AddRight(BitOperation::bitToIndex[current.OutEdges()]);

    next.SetNode(&nodes[hash_graph->GetNode(kmer)->Data()]);
    if (next.GetBeginKmer() != kmer)
        next.ReverseComplement();

    return next.InDegree() == 1;
}

