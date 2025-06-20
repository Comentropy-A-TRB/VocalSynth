// latest update on 2025/6/8

#include "Setting.h"

#include<cstdio>
#include<cstring>
#include<string>
#include<vector>
#include<iostream>
#include<io.h>
#include<tuple>

#include "AudioFile.h"
#include "Note.h"
#include "Singer.h"
#include "SoundTrans.h"
#include "Generator.h"

using std::ifstream,std::ofstream;
using std::vector,std::map,std::string;
using std::cin,std::cout;
// const int numChannels=1;
// const double sampleRate=44100;
// const float frequency=440.f;
// g++ .\SoundTrans.cpp main.cpp -o main -O2 -Wall -std=c++17 -L"./fftw/" -lfftw3-3 -lfftw3f-3 -lfftw3l-3
Singer singer;
struct Sheet{
    vector<Note> Que;
}Tar;
string pathOfSheet,dictPath;

void getBasicInfo(){
    std::cerr<<"please input dictionary path first. And then BPM, path of sheet"<<std::endl;
    cin>>dictPath;  // like ./English_Dictionary/
    cin>>BPM>>pathOfSheet;
}
void getDictionary(){
    // Path of dictionary.
    _finddata_t file;
    intptr_t readRet;   // fixed! 之前我用 GCC8.1 的时候直接用 int 存了，现在改了 64 位指针标准，必须换类型，不然 RE。
    // std::cerr<<(dictPath+"*.txt")<<std::endl;
    if((readRet=_findfirst((dictPath+"*.txt").c_str(),&file))==-1)
        cout<<"Dictionary of singer is empty"<<std::endl,exit(0);
    string dic;int typ;
    // !!!!!!!!!!!!!!!!!!! 这里需要相当大的调整，具体来说是要把响度曲线直接加入可变，共振峰机制不太改变。
    // 另外关于共振峰拟合的问题，我们需要进行一番修正（在 singer 里面）。
    // 最后，我们可能还需要解决 —— 各种响度的适配问题。
    vector<int> F(4),T,loud;
    vector<double> Floud(4),Mix(3);
    do{
        if(DEBUG)cout<<"Find: "<<file.name<<std::endl;
        dic=file.name;
        ifstream pathofdic(dictPath+dic);
        assert(dic.size()>4);
        assert(dic.substr(dic.size()-4)==".txt");
        dic.resize(dic.size()-4);
        // 
        pathofdic>>typ;
        singer.setType(dic,typ);
        if(typ==3){
            string v1,v2;
            pathofdic>>v1>>v2;
            for(int i=0;i<3;i++)
                pathofdic>>Mix[i];
            singer.setMixVow(dic,v1,v2,Mix);
        }else{
            for(int i=1;i<=3;i++)
                pathofdic>>F[i]>>Floud[i];  // Note: change to 'Floud' !
            if(DEBUG)cout<<dic<<" "<<F[1]<<" "<<F[2]<<" "<<F[3]<<std::endl;
            singer.setFormant(dic,F,Floud);
        }
        int sz;
        pathofdic>>sz;T.resize(sz),loud.resize(sz);
        for(int i=0;i<sz;i++)    // ADSR is begin from zero
            pathofdic>>T[i]>>loud[i];
        assert(T.front()==0 && T.back()==100);
        singer.setADSR(dic,T,loud);
        //
    }while(_findnext(readRet,&file)==0);
    std::cerr<<"Dictionary got."<<std::endl;
}

template<typename T1,typename T2>
inline void readParameter(vector<pair<T1,T2> > &in,ifstream &readpath){
    int sz;readpath>>sz;in.resize(sz);
    for(auto &[x,y]:in)
        readpath>>x>>y;
}
void getSheetAndGen(){
    AudioFile<double>::AudioBuffer buff;
    ifstream readpath(pathOfSheet);
    string lyr,key;
    double beat,sumLength=0.0;
    vector<pair<double,double> > pit;
    vector<pair<int,int> > dyn,EQ;
    while(readpath>>key>>beat>>lyr){ 
        readParameter<double,double>(pit,readpath);
        readParameter<int,int>(dyn,readpath);
        readParameter<int,int>(EQ,readpath);
        Tar.Que.emplace_back((Note){lyr,key,beat,pit,dyn,EQ});
        sumLength+=beat;
    }
    buff.resize(numChannels);
    int totSamples=beatToSamples(sumLength,BPM,sampleRateGlobal)+2*vowDelaySamples;
    buff[0].resize(totSamples);
    int pos=0;
    for(Note &i:Tar.Que){
        int typ=singer.getType(i.lyr);
        if(typ==1)
            GenVow(pos,i,buff,true,false,singer);
        else if(typ==0)
            GenCon(pos,i,buff,false,singer);
        else if(typ==2){
            GenVow(pos,i,buff,false,true,singer);
            GenCon(pos,i,buff,true,singer);
        }else{
            GenMix(pos,i,buff,singer);
        }
    }
    AudioFile<double> audioFile;
    audioFile.setNumChannels(numChannels);
    audioFile.setSampleRate(sampleRateGlobal);
    audioFile.setAudioBuffer(buff);
    audioFile.save("Export.wav",AudioFileFormat::Wave);
}
int main(){
    freopen("DEBUG.txt","w",stdout);
    srand(time(NULL));
    getBasicInfo(); // get BPM, Path of Sheet, (channels, Bitdep, sampleRate.)
    getDictionary();    // get config of the singer.
    getSheetAndGen();
    // system("Pause");
    return 0;
} 

/*
./Japanese_Dictionary/
183 reyichang.txt
*/