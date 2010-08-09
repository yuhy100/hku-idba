#include "globals.h"
#include "Sequence.h"
#include "Log.h"
#include "Reader.h"
#include "Writer.h"
#include "Utils.h"

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <cstring>
#include <set>
#include <string>
#include <map>

using namespace std;

const int MaxSeqLength = 10000000;

int reverse_table[256] = {0};
char buf[MaxLine];
char line[MaxLine];
char readName[MaxLine];
char contigName[MaxLine];
char comment[MaxLine];
char strand[MaxLine];
vector<string> names;
vector<Sequence> contigs;
map<string, int> dict;
vector<vector<int> > matchs;
long long errors[100] = {0};
long long partial[100] = {0};
int cover_count[10000] = {0};

int main(int argc, char *argv[])
{
    reverse_table['A'] = 'T';
    reverse_table['C'] = 'G';
    reverse_table['G'] = 'C';
    reverse_table['T'] = 'A';
    reverse_table['N'] = 'N';

    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: calcCover contig-file read-file blat-file\n");
        fprintf(stderr, "       [--minContig l]\n");
        exit(1);
    }

    FastAReader contigs_reader(argv[1]);
    Sequence seq;
    string comment;
    while (contigs_reader.Read(seq, comment))
    {
        sscanf(comment.c_str(), "%s", contigName);
        dict[contigName] = contigs.size();
        contigs.push_back(seq);
        names.push_back(contigName);
    }

    matchs.resize(contigs.size());
    for (unsigned i = 0; i < contigs.size(); ++i)
        matchs[i].resize(contigs[i].Size(), 0);

    FastAReader reader(argv[2]);
    FILE *fp = fopen(argv[3], "rb");

    fgets(line, MaxLine, fp);
    fgets(line, MaxLine, fp);
    fgets(line, MaxLine, fp);
    fgets(line, MaxLine, fp);
    fgets(line, MaxLine, fp);

    Sequence current_read;
    string current_read_name;
    long long total_aligned_reads = 0;

    while (fgets(line, MaxLine, fp) != NULL)
    {
        int match_count = 0;
        int read_from, read_to;
        int contig_from, contig_to;
        int contig_length;
        sscanf(line, "%d %s %s %s %s %s %s %s %s %s %s %d %d %s %d %d %d",
                &match_count, buf, buf, buf, buf, buf, buf, buf, strand,
                readName, buf, &read_from, &read_to, contigName, &contig_length, &contig_from, &contig_to
              );

        if (current_read_name == readName)
            continue;

        ++total_aligned_reads;

        while(current_read_name != readName)
            reader.Read(current_read, current_read_name);

        read_from;
        read_to--;
        contig_from;
        contig_to--;

        //if (read_from > read_to)
//        if (strand[0] == '-')
//        {
////            int from = current_read.Size() - 1 - read_to;
////            int to = current_read.Size() - 1 - read_from;
////            read_from = from;
////            read_to = to;
//
////            int from, to;
////            from = contig_length - 1 - contig_to;
////            to = contig_length - 1 - contig_from;
////            contig_from = from;
////            contig_to = to;
//        }

        int index = dict[contigName];

        int step = 1;
//        if (contig_from > contig_to)
//            step = -1;
        if (strand[0] == '-')
            step = -1;

        bool is_prev_error = false;
        //cout << current_read << endl;
        for (int i = 0; i < current_read.Size(); ++i)
        {
//            if (i < read_from || i > read_to)
//                ++errors[i];
//            else
            {
                int offset = i - read_from;
                int contig_index = contig_from + step * offset;

                if (strand[0] == '-')
                    contig_index = contig_to + step * offset;
                
                //cout << contigs[index][contig_index] << " " << contig_index << " ";
                if (strand[0] == '+')
                {
                    if (contigs[index][contig_index] != current_read[i])
                    {
                        ++errors[i];
                        if (is_prev_error)
                            ++partial[i];
                        is_prev_error = true;
                    }
                    else
                        is_prev_error = false;
                }
                else
                {
                    if (contigs[index][contig_index] != reverse_table[current_read[i]])
                    {
                        ++errors[i];
                        if (is_prev_error)
                            ++partial[i];
                        is_prev_error = true;
                    }
                    else
                        is_prev_error = false;
                }
            }
        }
        //cout << endl << endl;

        if (contig_from > contig_to)
            swap(contig_from, contig_to);

        for (int i = contig_from; i <= contig_to; ++i)
            matchs[index][i]++;
    }

    long long sum_matchs = 0;
    long long sum_lengths = 0;
    for (unsigned i = 0; i < matchs.size(); ++i)
    {
        sum_lengths += matchs[i].size();

        long long sum = 0;
        for (unsigned j = 0; j < matchs[i].size(); ++j)
        {
            sum += matchs[i][j];
            //++cover_count[matchs[i][j]];
        }

        sum_matchs += sum;

        ++cover_count[sum/matchs[i].size()];
    }
    printf("coverage %.2f\n", 1.0 * sum_matchs / sum_lengths);

    int total_reads = reader.NumReads();
    cout << total_aligned_reads << endl;
    for (int i = 1; i < 100; ++i)
        cout << i << " " << errors[i] << " " << errors[i] * 1.0 / total_aligned_reads << " " << partial[i] * 1.0 / errors[i-1] << endl;

    for (int i = 0; i < 200; ++i)
        cout << i << " " << cover_count[i] << endl;

    return 0;
}

