#include "globals.h"
#include "AbstractAssembler.h"

#include <vector>

using namespace std;


AbstractAssembler::AbstractAssembler(const Option &option, int64 read_buffer)
{
    this->option = option; 
    this->read_buffer = read_buffer;
    read_num = 0;
    read_reader = new FastAReader(option.readfile);
}

vector<Read> &AbstractAssembler::GetShortReads()
{
    if (read_buffer == 0)
    {
        if (read_num != 0)
            return empty_reads;

        if (reads.size() == 0)
            ReadBunchOfReads(reads, read_reader->NumReads());
        
        read_num = reads.size();
        return reads;
    }
    else
    {
        ReadBunchOfReads(reads, read_buffer);
        read_num += reads.size();
        return reads;
    }
}

void AbstractAssembler::ReadBunchOfReads(vector<Read> &reads, int64 expected_num)
{
    reads.resize(expected_num);
    Sequence seq;
    string comment;
    int count = 0;
    for (int64 i = 0; i < (int64)reads.size(); ++i)
    {
        if (!read_reader->Read(seq, comment))
        {
            reads.resize(i);
            break;
        }

        reads[i].Inactivate();
        seq.Trim(option.trim);

        int from = 0;
        int to = 0;
        int last = 0;
        
        for (int current = 0; current < seq.Size(); ++current)
        {
            if (seq[current] == 'N')
            {
                if (current - last > to - from)
                {
                    from = last;
                    to = current;
                }

                last = current + 1;
            }
        }

        if (seq.Size() - last > to - from)
        {
            from = last;
            to = seq.Size();
        }

        Sequence valid_seq;
        seq.GetSubSequence(valid_seq, from, to - from);

        seq = valid_seq;

//        seq.TrimError();
//        seq.ReverseComplement();
//        seq.TrimError();
//        seq.ReverseComplement();


        if (!seq.IsChar())
            continue;

        ++count;

        seq.Encode();
        reads[i] = seq;
    }

    printf("reads %d\n", count);
}

