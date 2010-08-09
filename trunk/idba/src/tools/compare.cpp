#include <iostream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

vector<int> record1;
vector<int> record2;
vector<int> lengths1;
vector<int> lengths2;

int matchPair(int offset1, int length1, int offset2, int length2)
{
    if (offset1 == -1 || offset2 == -1)
        return -1;

    if (length1 < 100 || length2 < 100)
        return -1;

    if (offset1 < offset2)
    {
        if (abs(offset1 + length1 - offset2) < 500)
            return 0;
        else
            return -1;
    }
    else
    {
        if (abs(offset2 + length2 - offset1) < 500)
            return 1;
        else
            return -1;
    }
}

int main(int argc, char *argv[])
{
    ifstream fin1(argv[1], ios_base::in | ios_base::binary);
    ifstream fin2(argv[2], ios_base::in | ios_base::binary);

    string line;
    fin1 >> line;

    int size;
    fin1 >> size;

    fin2 >> line;
    fin2 >> line;

    record1.resize(size);
    lengths1.resize(size);
    record2.resize(size);
    lengths2.resize(size);

    int id, offset, length;
    while (fin1 >> id >> offset >> length)
    {
        //cout << id << " " << offset << endl;
        record1[id] = offset;
        lengths1[id] = length;
    }

    while (fin2 >> id >> offset >> length)
    {
        record2[id] = offset;
        lengths2[id] = length;
    }

    //cout << size << endl;

    int count1 = 0;
    int count2 = 0;
    int match = 0;
    for (int i = 0; i < size; ++i)
    {
        for (int j = i+1; j < size; ++j)
        {
            //if (record1[i] != -1 && record1[j] != -1)
            if (matchPair(record1[i], lengths1[i], record1[j], lengths1[j]) != -1)
                ++count1;

            //if (record2[i] != -1 && record2[j] != -1)
            if (matchPair(record2[i], lengths2[i], record2[j], lengths2[j]) != -1)
                ++count2;

            if (matchPair(record1[i], lengths1[i], record1[j], lengths1[j]) != -1 &&
                matchPair(record2[i], lengths2[i], record2[j], lengths2[j]) != -1)
            {
                if (matchPair(record1[i], lengths1[i], record1[j], lengths1[j]) == matchPair(record2[i], lengths2[i], record2[j], lengths2[j]))
                    ++match;
            }
        }
    }

    printf("%d,%d,%d,%.6f,%.6f\n", count1, count2, match, match*1.0/count1, match*1.0/count2);

    return 0;
}
