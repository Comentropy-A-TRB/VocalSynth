#define _USE_MATH_DEFINES
#include<cstdio>
#include<cstring>
#include<string>
#include<vector>
#include<cmath>
#include<map>
#include<iostream>
#include<random>
#include<ctime>
#include<io.h>
#include<tuple>
#include "AudioFile.h"
#include "SoundTrans.h"
#include "Singer.h"
#define DEBUG 1
#define vowDelaySamples 50
#define conDelaySamples 10
#define eps 1e-6
using std::ifstream,std::ofstream;
using std::vector,std::map,std::string;
using std::cin,std::cout;
const int numChannels=1;
const double sampleRate=44100;
// const float frequency=440.f;
// g++ .\SoundTrans.cpp main.cpp -o main -O2 -Wall -std=c++17 -L"./fftw/" -lfftw3-3 -lfftw3f-3 -lfftw3l-3
Singer singer;
struct Note{
    string lyr,key;
    double beat;
};
struct Sheet{
    vector<Note> Que;
}Tar;
int BPM;
string pathOfSheet;
void getBasicInfo(){
    cin>>BPM>>pathOfSheet;
}
void getDictionary(){
    // Path of dictionary.
    _finddata_t file;
    int readRet;
    if((readRet=_findfirst("./Dictionary/*.txt",&file))==-1)
        cout<<"Dictionary of singer is empty"<<std::endl,exit(0);
    string dic;int typ;
    vector<int> F(4),T(4),EQ(4);
    vector<double> loud(4),Mix(3);
    do{
        if(DEBUG)cout<<"Find: "<<file.name<<std::endl;
        dic=file.name;
        ifstream pathofdic("./Dictionary/"+dic);
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
                pathofdic>>F[i]>>loud[i];
            if(DEBUG)cout<<dic<<" "<<F[1]<<" "<<F[2]<<" "<<F[3]<<std::endl;
            singer.setFormant(dic,F,loud);
        }
        for(int i=0;i<4;i++)    // ADSR is begin from zero
            pathofdic>>T[i]>>EQ[i];
        assert(T[0]==0 && T[3]==100);
        singer.setADSR(dic,T,EQ);
        //
    }while(_findnext(readRet,&file)==0);
}
void realizeADSR(int pos,int totsamp,Note note,vector<double> &buff){
    for(int i=0;i<totsamp;i++){
        buff[pos+i]*=singer.getADSR(note.lyr,floor(100.0*i/totsamp));
    }
}
void GenCon(int &pos,Note note,AudioFile<double>::AudioBuffer &buff,bool KeepFormant=false){
    if(DEBUG)cout<<"Generating note: "<<note.lyr<<std::endl;
    int totsamp=beatToSamples(note.beat,BPM,sampleRate)+conDelaySamples;
    // 湍流就是白噪音加上均衡。
    const double mean = 0.0;//均值
    const double stddev = 0.1;//标准差
    std::mt19937 rd(rand());
    std::normal_distribution<double> dist(mean, stddev);
    for(int i=0;i<totsamp;i++){
        int to=pos+i;
        for(int channel=0;channel<numChannels;channel++){
            buff[channel][to]+=std::min(0.95,std::max(-0.95,dist(rd)));
        }
    }
    vector<std::tuple<double,double,double> > req;
    if(!KeepFormant)
        req.push_back({0,2000,0});
    for(int i=2000;i<180000;i+=20){
        req.push_back({i,i+20,std::max(.0,singer.getConCoeffic(note.lyr,i)/3)});
    }
    req.push_back({18000,20000,0});
    for(int channel=0;channel<numChannels;channel++)
        FFTfilter(pos,totsamp,buff[channel],buff[channel],sampleRate,req,1);
    for(int channel=0;channel<numChannels;channel++)
        realizeADSR(pos,totsamp,note,buff[channel]);
    pos+=totsamp-conDelaySamples;
}
void GenVow(int &pos,Note note,AudioFile<double>::AudioBuffer &buff,bool updPos=true,bool NoFormant=false){
    if(DEBUG)cout<<"Generating note: "<<note.lyr<<std::endl;
    int totsamp=beatToSamples(note.beat,BPM,sampleRate)+vowDelaySamples;  // 合成后置衔接过渡
    auto genAccordingFreq=[&](double f0,double k){
        if(NoFormant)
            k*=0.08;
        else
            k*=std::max(eps,singer.getCoefficient(note.lyr,f0));
        if(DEBUG)cout<<"req: "<<f0<<" "<<k<<std::endl;
        for(int i=0;i<totsamp;i++){
            int to=pos+i;
            for(int channel=0;channel<numChannels;channel++){
                buff[channel][to]+=k*sin((static_cast<double>(i)/sampleRate)*f0*2.0*M_PI);
            }
        }
    };
    double freq=keyToFreq(note.key);
    double now=freq;
    if(note.lyr=="SP")  goto finished;
    for(int i=0;i<20;i++){ 
        // 生成 10 条 泛音，从第三条之后大小线性递减。
        // upd: 由于已经规定了共振峰和二次函数拟合，所以可以不用递减了。
        // if(i>3)
        //     k-=0.04;
        genAccordingFreq(now,1.0/20);    // 记得平均。
        now+=freq;
    }
    for(int channel=0;channel<numChannels;channel++)
        realizeADSR(pos,totsamp,note,buff[channel]);
    finished:
    if(DEBUG)cout<<"Generated!"<<std::endl;
    if(updPos)
        pos+=totsamp-vowDelaySamples;
}
void GenMix(int &pos,Note note,AudioFile<double>::AudioBuffer &buff){
    if(DEBUG)cout<<"Generating note: "<<note.lyr<<std::endl;
    int totsamp=beatToSamples(note.beat,BPM,sampleRate)+vowDelaySamples;  // 合成后置衔接过渡
    auto &X=singer.Mix[note.lyr];
    auto genAccordingFreq=[&](double f0,double k){
        double k1=singer.getCoefficient(X.v1,f0),k2=singer.getCoefficient(X.v2,f0);
        double beg=X.tm1*(1-X.mixK),ed=X.tm2+X.tm1+X.tm2*X.mixK,nowK;
        beg*=totsamp,ed*=totsamp;
        for(int i=0;i<totsamp;i++){
            int to=pos+i;
            if(i<beg)   nowK=k1;
            else if(i<ed)   nowK=(k2-k1)/(ed-beg)*(i-beg)+k1;
            else    nowK=k2;
            for(int channel=0;channel<numChannels;channel++){
                buff[channel][to]+=k*nowK*sin((static_cast<double>(i)/sampleRate)*f0*2.0*M_PI);
            }
        }
    };
    double freq=keyToFreq(note.key);
    double now=freq;
    for(int i=0;i<20;i++){ 
        genAccordingFreq(now,1.0/20);    // 记得平均。
        now+=freq;
    }
    for(int channel=0;channel<numChannels;channel++)
        realizeADSR(pos,totsamp,note,buff[channel]);
    if(DEBUG)cout<<"Generated!"<<std::endl;
    pos+=totsamp-vowDelaySamples;
}
void getSheetAndGen(){
    AudioFile<double>::AudioBuffer buff;
    ifstream readpath(pathOfSheet);
    string lyr,key;double beat,sumLength=0.0;
    while(readpath>>key>>beat>>lyr){ 
        Tar.Que.push_back((Note){lyr,key,beat});
        sumLength+=beat;
    }
    buff.resize(numChannels);
    int totSamples=beatToSamples(sumLength,BPM,sampleRate)+2*vowDelaySamples;
    buff[0].resize(totSamples);
    int pos=0;
    for(Note &i:Tar.Que){
        int typ=singer.getType(i.lyr);
        if(typ==1)
            GenVow(pos,i,buff);
        else if(typ==0)
            GenCon(pos,i,buff);
        else if(typ==2){
            GenVow(pos,i,buff,false,true);
            GenCon(pos,i,buff,true);
        }else{
            GenMix(pos,i,buff);
        }
    }
    AudioFile<double> audioFile;
    audioFile.setNumChannels(1);
    audioFile.setSampleRate(sampleRate);
    audioFile.setAudioBuffer(buff);
    audioFile.save("Export.wav",AudioFileFormat::Wave);
}
int main(){
    srand(time(NULL));
    getBasicInfo(); // get BPM, Path of Sheet, (channels, Bitdep, sampleRate.)
    getDictionary();    // get config of the singer.
    getSheetAndGen();
    // system("Pause");
    return 0;
}