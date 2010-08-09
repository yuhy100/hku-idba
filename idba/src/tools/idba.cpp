#include "config.h"

#include "globals.h"
#include "Log.h"
#include "Sequence.h"
#include "SimpleAssembler.h"
#include "IterativeAssembler.h"
#include "ScaffoldAssembler.h"
#include "GappedScaffoldAssembler.h"
#include "Option.h"
#include "Reader.h"
#include "Writer.h"

#include <unistd.h>
#include <sys/wait.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace std;

Option option;

void BuildSimpleGraph()
{
    SimpleAssembler *simple = new SimpleAssembler(option);
    vector<Contig> *contigs = new vector<Contig>();
    simple->Assemble(*contigs);
    delete simple;

    FastAWriter writer(option.graphfile);
    for (unsigned i = 0; i < (*contigs).size(); ++i)
    {
        (*contigs)[i].Decode();
        writer.WriteFormat((*contigs)[i], "contig%d", i);
    }

    delete contigs;
}

void Iterate()
{
    IterativeAssembler simple(option);
    vector<Contig> contigs;
    simple.Assemble(contigs);

    FastAWriter writer(option.contigfile);
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        contigs[i].Decode();
        writer.WriteFormat(contigs[i], "contig%d", i);
    }
}

void Scaffold()
{
    kmerLength = option.maxk;
    ScaffoldAssembler simple(option);
    vector<Contig> contigs;
    simple.Assemble(contigs);

    FastAWriter writer(option.scaffile);
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        contigs[i].Decode();
        writer.WriteFormat(contigs[i], "contig%d", i);
    }
}

void ScaffoldWithGap()
{
    GappedScaffoldAssembler simple(option);
    vector<Contig> contigs;
    simple.Assemble(contigs);

    FastAWriter writer(option.gap_scaffold_file);
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        contigs[i].Decode();
        writer.WriteFormat(contigs[i], "contig%d", i);
    }
}

void Usage()
{
    cout << "IDBA: Iterative De Bruijn graph short read Assembler" << endl;
    cout << "Version " << VERSION << endl;
    cout << endl;
    cout << "Usage: idba --read read-file [--output out] [options]\n" << endl;
}

int main(int argc, char *argv[])
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("read,r", po::value<string>(&option.readfile), "read file")
        ("read2", po::value<string>(&option.readfile2), "read file")
        ("output,o", po::value<string>(&option.prefix)->default_value("out"), "prefix of output")
        ("scaffold", "use pair end information to merge contigs")
        ("gap", "use pair end information to merge contigs with gap")
        ("mink", po::value<int>(&option.mink)->default_value(25), "minimum k value")
        ("maxk", po::value<int>(&option.maxk)->default_value(50), "maximum k value")
        ("minCount", po::value<int>(&option.min_count)->default_value(2), "filtering threshold for each k-mer")
        ("cover", po::value<double>(&option.cover)->default_value(0), "the cutting coverage for contigs")
        ("trim", po::value<int>(&option.trim)->default_value(0), "discard n bases a the end of each read")
        ("minPairs", po::value<int>(&option.min_pairs)->default_value(5), "minimum number of pair-end connections to join two contigs")
        ("minLength", po::value<int>(&option.min_length)->default_value(100), "minimum contig length for output")
    ;

    try
    {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        option.is_scaffold = vm.count("scaffold");
        option.is_gap = vm.count("gap");

        if (vm.count("help") || vm.count("read") == 0)
        {
            Usage();
            cout << desc << endl;
            return 1;
        }
    }
    catch (exception &e)
    {
        Usage();
        cout << desc << endl;
        return 1;
    }

    option.Compute();

    if (!fs::exists(option.graphfile))
    {
        int pid = fork();
        if (pid == 0)
        {
            BuildSimpleGraph();
            exit(1); 
        }
        else
        {
            int status;
            waitpid(pid, &status, 0);
        }
    }

    if (!fs::exists(option.contigfile))
    {
        int pid = fork();
        if (pid == 0)
        {
            Iterate();
            exit(1); 
        }
        else
        {
            int status;
            waitpid(pid, &status, 0);
        }
    }

    if (option.is_scaffold && !fs::exists(option.scaffile))
    {
        int pid = fork();
        if (pid == 0)
        {
            Scaffold();
            exit(1); 
        }
        else
        {
            int status;
            waitpid(pid, &status, 0);
        }
    }

    if (option.is_scaffold && option.is_gap && !fs::exists(option.gap_scaffold_file))
    {
        int pid = fork();
        if (pid == 0)
        {
            ScaffoldWithGap();
            exit(1); 
        }
        else
        {
            int status;
            waitpid(pid, &status, 0);
        }
    }

    return 0;
}

