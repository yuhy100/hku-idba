#include "globals.h"
#include "Sequence.h"
#include "Utils.h"
#include "Reader.h"

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <set>
#include <string>

using namespace std;

char buf[MaxLine];
char line[MaxLine];
char name[MaxLine];
vector<int> recrod;
vector<int> lengths;
vector<int> valid_lengths;
set<string> valid_contigs;
int min_contig = 100;
double similar = 0.99;

int main(int argc, char *argv[])
{
    AddParameter("minContig", &min_contig, INTEGER);
    AddParameter("similar", &similar, FLOAT);

    ProcessParameters(argc, argv);

    //cout << min_contig << " " << similar << endl;

    FastAReader ref_reader(argv[1]);
    FastAReader qry_reader(argv[2]);

    Sequence ref;
    string comment;
    ref_reader.Read(ref, comment);

    Sequence seq;
    string tmp;
    int index = 0;
    while (qry_reader.Read(seq, tmp))
        ++index;
    recrod.resize(index);
    lengths.resize(index);
    fill(recrod.begin(), recrod.end(), -1);
    fill(lengths.begin(), lengths.end(), -1);

    vector<bool> flags(ref.Size() + 1);

    while (fgets(line, MaxLine, stdin) != NULL)
    {
        if (line[0] == '>')
        {
            int ref_len;
            int contig_len;
            sscanf(line, "%s %s %d %d", buf, name, &ref_len, &contig_len);

            int id = atoi(name + 6);

            while (fgets(line, MaxLine, stdin) != NULL)
            {
                if (line[0] == '>')
                {
                    fseek(stdin, -strlen(line), SEEK_CUR);
                    break;
                }

                int ref_from, ref_to;
                int contig_from, contig_to;
                int errors;
                sscanf(line, "%d %d %d %d %d", &ref_from, &ref_to, &contig_from, &contig_to, &errors);

                --ref_from;
                --ref_to;
                --contig_from;
                --contig_to;

                if (ref_from > ref_to)
                    swap(ref_from, ref_to);
                ++ref_to;

                if (contig_from > contig_to)
                    swap(contig_from, contig_to);
                ++contig_to;

                int ref_match_len = ref_to - ref_from;
                int contig_match_len = contig_to - contig_from;

                if (ref_match_len >= similar * contig_len
                        && contig_len >= ref_match_len * similar
                        && errors <= contig_len * (1 - similar)
                        && contig_len >= min_contig)
                //if (contig_len >= 100)
                {
                    for (int i = ref_from; i < ref_to; ++i)
                        flags[i] = 1;

                    //cout << id << " " << ref_from << endl;
                    recrod[id] = ref_from;
                    lengths[id] = contig_match_len;
                    valid_lengths.push_back(contig_match_len);
                    valid_contigs.insert(name);
                }

                while (fgets(line, MaxLine, stdin) != NULL)
                {
                    if (line[0] == '0')
                        break;
                }
            }
        }
    }

    int count = 0;
    for (int i = 0; i < flags.size(); ++i)
    {
        if (flags[i])
            ++count;
    }

    sort(valid_lengths.begin(), valid_lengths.end());
    reverse(valid_lengths.begin(), valid_lengths.end());
    int n50 = 0;
    int sum = 0;


    for (unsigned i = 0; i < valid_lengths.size(); ++i)
    {
        sum += valid_lengths[i];
        if (sum >= 0.5 * ref.Size() && n50 == 0)
            n50 = valid_lengths[i];
    }

    int maximum = 0;
    int mean = 0;
    if (valid_lengths.size() > 0)
    {
        maximum = valid_lengths[0];
        mean = sum / valid_lengths.size();
    }

    qry_reader.Rewind();
    int num_contigs = 0;
    long long sum_wrong = 0;
    int num_wrong = 0;
    while (qry_reader.Read(seq, comment))
    {
        if (seq.Size() >= min_contig)
        {
            ++num_contigs;

            sscanf(comment.c_str(), "%s", name);
            if (valid_contigs.find(name) == valid_contigs.end())
            {
                ++num_wrong;
                sum_wrong += seq.Size();
//                cout << ">" << comment << endl;
//                cout << seq << endl;
            }
            else
            {
            }
        }
    }

    printf("contigs: %d N50: %d coverage: %.2f%% max: %d mean:%d\n", 
            (int)valid_contigs.size(), n50, count * 100.0 / flags.size(),
            maximum, mean);
    printf("error: %d %lld\n", num_wrong, sum_wrong);

    return 0;
}
