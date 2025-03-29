// latest update on 2025/3/28
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
#define vowDelaySamples 25
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
    vector<pair<double,double> > pit;   // 内部保存采用比值保存，相对音高记法。
    vector<pair<int,int> > dyn; // 采用比值保存，相对响度记法（百分数），如果 !empty 则替换默认 ADSR。
    vector<pair<int,int> > EQ;  // 直接对频域保存，限制频域为 20k，整数下标、相对响度（百分数）记法。
};
struct Sheet{
    vector<Note> Que;
}Tar;
int BPM;
string pathOfSheet;
namespace MathCalc{
    double getLinearVal(double x0,double y0,double x1,double y1,double x){
        if(fabs(x0-x1)<eps) return (y0+y1)/2;
        return 1.0*(y0-y1)/(x0-x1)*(x-x0)+y0;
    }
    double Mod(double x,double y){
        return x-floor(x/y)*y;  // 笑死我了没写 floor.
    }
}
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
    // !!!!!!!!!!!!!!!!!!! 这里需要相当大的调整，具体来说是要把响度曲线直接加入可变，共振峰机制不太改变。
    // 另外关于共振峰拟合的问题，我们需要进行一番修正（在 singer 里面）。
    // 另外我们可能需要写一个 Math 库了（现在不写，还够用），来解决直线拟合/曲线拟合的问题，方便计算。
    // 最后，我们可能还需要解决 —— 各种响度的适配问题。
    vector<int> F(4),T,loud;
    vector<double> Floud(4),Mix(3);
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
}
void realizeDefaultADSR(int pos,int totsamp,const Note &note,vector<double> &buff){
    for(int i=0;i<totsamp;i++){
        buff[pos+i]*=singer.getADSR(note.lyr,floor(100.0*i/totsamp));
    }
    // 注意先做这个 ADSR 再调响度，为了比较方便调教，把元音的默认 ADSR 调得保守一点，形状如一个梯形
}
void updateDyn(int pos,int totsamp,const Note &note,vector<double> &buff){
    auto &dyn=note.dyn;
    if(dyn.empty()) return ;
    for(int i=0,j=0,sz=dyn.size();i<totsamp;i++){
        double perc=100.0*i/totsamp;
        while(j+1<sz && perc>=dyn[j+1].first)
            ++j;
        buff[pos+i]*=MathCalc::getLinearVal(dyn[j].first,dyn[j].second,dyn[j+1].first,dyn[j+1].second,perc)/100.0;  // 注意这里的百分运算。
    }
}
void updateEQ(int pos,int totsamp,const Note &note,vector<double> &buff){
    auto &EQ=note.EQ;
    if(EQ.empty())  return ;
    vector<std::tuple<double,double,double> > req;
    for(int i=0,j=0,sz=EQ.size();i<20000;i+=10){
        while(j+1<sz && i>=EQ[j+1].first)
            ++j;
        req.push_back({i,i+10,MathCalc::getLinearVal(EQ[j].first,EQ[j].second,EQ[j+1].first,EQ[j+1].second,i)*0.01});
    }
    FFTfilter(pos,totsamp,buff,buff,sampleRate,req,1);
}
void setDynZero(int pos,int totsamp,vector<double> &buff){
    // 考虑音量均一化算法。（必定要保证波形和原来差不多地把音量提升）
}
void effectParameters(int pos,int totsamp,const Note &note,AudioFile<double>::AudioBuffer &buff){
    // 注意因为各种音的 Delay 不同，所以这里把 totsamp 传进来。
    // 提前做一下响度均一化。
    for(int channel=0;channel<numChannels;channel++){
        setDynZero(pos,totsamp,buff[channel]);
        realizeDefaultADSR(pos,totsamp,note,buff[channel]);
        updateDyn(pos,totsamp,note,buff[channel]);
        updateEQ(pos,totsamp,note,buff[channel]);
    }
}
void GenCon(int &pos,const Note &note,AudioFile<double>::AudioBuffer &buff,bool KeepFormant=false){
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
    for(int i=2000;i<18000;i+=20){  // 18k????? fixed.
        req.push_back({i,i+20,std::max(.0,singer.getConCoeffic(note.lyr,i)/3)});
    }
    req.push_back({18000,20000,0});
    for(int channel=0;channel<numChannels;channel++)
        FFTfilter(pos,totsamp,buff[channel],buff[channel],sampleRate,req,1);
   
    effectParameters(pos,totsamp,note,buff);
        
    // for(int channel=0;channel<numChannels;channel++)
    //     realizeDefaultADSR(pos,totsamp,note,buff[channel]);
    pos+=totsamp-conDelaySamples;
}
void GenVow(int &pos,const Note &note,AudioFile<double>::AudioBuffer &buff,bool updPos=true,bool NoFormant=false){
    if(DEBUG)cout<<"Generating note: "<<note.lyr<<std::endl;
    int totsamp=beatToSamples(note.beat,BPM,sampleRate)+vowDelaySamples;  // 合成后置衔接过渡
    auto &pit=note.pit;
    auto genAccordingFreq=[&](double f0,double k,int lay){  // upd: lay - layer 表示音高偏差倍数
        if(NoFormant)  k*=0.08;
        if(DEBUG)cout<<"req: "<<f0<<" "<<k<<" "<<pit.size()<<std::endl;
        double w=(pit.empty()?f0:pit[0].second*lay+f0)*2.0*M_PI,b=0;   // 第一个位置就是零号点，而且此时相位是 0. 整的
        // 不要忘了 *lay+f0
        for(int i=0,j=0,pitsz=pit.size();i<totsamp;i++){
            int to=pos+i;
            while(j+1<pitsz && i>=pit[j+1].first*totsamp)
                ++j;
            if(pitsz && i>0){
                double now=MathCalc::getLinearVal(pit[j].first*totsamp,pit[j].second,   // 注意不要给 second*totsamp
                                         pit[j+1].first*totsamp,pit[j+1].second,i)*lay+f0;   // 音高偏差加 f0 得到真实音高
                // if(DEBUG && i%100==0)  cout<<now<<" "<<w<<" "<<b<<" "<<i<<"/"<<totsamp<<std::endl;
                // 计算相位使得两者在上一个点重合。并保证它最简。（不会太大）
                // w*(i-1) + b = w'*(i-1)+ b'
                if(fabs(w-now*2.0*M_PI)>eps*40){   // 规避一些精度问题
                    b=MathCalc::Mod((w-now*2.0*M_PI)*(i-1)/sampleRate+b,2.0*M_PI);  // 我就说这里 b 怎么一直是 0，原来是上面没 floor
                    // Note：观察波形可以知道，后面这段完全不匀称，显然是相位有了问题。
                    w=now*2.0*M_PI;
                }
            }
            for(int channel=0;channel<numChannels;channel++){
                buff[channel][to]+=k*sin((static_cast<double>(i)/sampleRate)*w+b);
            }
            // NOTE：如果在这里强制使用之前的音高不变的合成方法就会爆炸，只能上 EQ 了。
        }
    };
    auto FormantEQ=[&](){   // 我就不信有人能发出 50 Hz 以下的声音。
        vector<std::tuple<double,double,double> > req;
        req.push_back({0,50,.0});
        for(int i=50;i<18000;i+=20)
            req.push_back({i,i+20,std::max(eps,singer.getCoefficient(note.lyr,i))});
        req.push_back({18000,20000,.0});
        for(int channel=0;channel<numChannels;channel++)
            FFTfilter(pos,totsamp,buff[channel],buff[channel],sampleRate,req,1);
    };
    double freq=keyToFreq(note.key);
    double now=freq;
    if(note.lyr=="SP")  goto finished;
    for(int i=0;i<20 && now<=18000;i++){ 
        // 生成 10 条 泛音，从第三条之后大小线性递减。
        // upd: 由于已经规定了共振峰和二次函数拟合，所以可以不用递减了。
        // if(i>3)
        //     k-=0.04;
        genAccordingFreq(now,1.0/20,i+1);    // 记得平均。
        now+=freq;  
    }
    FormantEQ();
    effectParameters(pos,totsamp,note,buff);
    // for(int channel=0;channel<numChannels;channel++)
    //     realizeDefaultADSR(pos,totsamp,note,buff[channel]);
    finished:
    if(DEBUG)cout<<"Generated!"<<std::endl;
    if(updPos)
        pos+=totsamp-vowDelaySamples;
}
void GenMix(int &pos,const Note &note,AudioFile<double>::AudioBuffer &buff){
    if(DEBUG)cout<<"Generating note: "<<note.lyr<<std::endl;
    int totsamp=beatToSamples(note.beat,BPM,sampleRate)+vowDelaySamples;  // 合成后置衔接过渡
    auto &X=singer.Mix[note.lyr];
    auto genAccordingFreq=[&](double f0,double k){
        double k1=singer.getCoefficient(X.v1,f0),k2=singer.getCoefficient(X.v2,f0);
        double beg=X.tm1*(1-X.mixK),ed=X.tm1+X.tm2*X.mixK,nowK;
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
    
    effectParameters(pos,totsamp,note,buff);
    // for(int channel=0;channel<numChannels;channel++)
    //     realizeDefaultADSR(pos,totsamp,note,buff[channel]);
    if(DEBUG)cout<<"Generated!"<<std::endl;
    pos+=totsamp-vowDelaySamples;
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
    freopen("DEBUG.txt","w",stdout);
    srand(time(NULL));
    getBasicInfo(); // get BPM, Path of Sheet, (channels, Bitdep, sampleRate.)
    getDictionary();    // get config of the singer.
    getSheetAndGen();
    // system("Pause");
    return 0;
} 