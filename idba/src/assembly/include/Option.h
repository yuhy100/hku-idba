#ifndef __OPTION_H_

#define __OPTION_H_

#include <string>

struct Option
{
    std::string prefix;
    std::string readfile;
    std::string readfile2;
    std::string graphfile;
    std::string contigfile;
    std::string scaffile;
    std::string gap_scaffold_file;
    int mink;
    int maxk;
    int min_length;
    unsigned prefix_length;
    int min_count;
    int trim;
    double cover;
    int min_pairs;
    int min_contig;
    bool is_scaffold;
    bool is_gap;
    uint64 mask;

    Option()
    {
        prefix = "out";
        mink = 25;
        maxk = 50;
        prefix_length = 3;
        min_count = 2;
        trim = 0;
        cover = 0;
        min_pairs = 5;
        is_scaffold = false;
        is_gap = false;
    }

    void Compute()
    {
        graphfile = prefix + ".graph";
        contigfile = prefix + "-contig.fa";
        scaffile = prefix + "-contig-mate.fa";
        gap_scaffold_file = prefix + "-contig-mate-gap.fa";

        mask = (1 << prefix_length) - 1;
    }
};

#endif
