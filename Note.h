#ifndef NOTE_H
#define NOTE_H 1

#include<vector>
#include<string>
#include<cstdio>

using std::vector;
using std::string;
using std::pair;

struct Note{
    string lyr,key;
    double beat;
    vector<pair<double,double> > pit;   // 内部保存采用比值保存，相对音高记法。
    vector<pair<int,int> > dyn; // 采用比值保存，相对响度记法（百分数），如果 !empty 则替换默认 ADSR。
    vector<pair<int,int> > EQ;  // 直接对频域保存，限制频域为 20k，整数下标、相对响度（百分数）记法。
};

#endif