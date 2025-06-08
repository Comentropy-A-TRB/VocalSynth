#ifndef GENERATOR_H
#define GENERATOR_H

#include "SoundTrans.h"
#include "Setting.h"
#include "Note.h"
#include "Effects.h"
#include "AudioFile.h"
#include "Singer.h"

#include<iostream>
#include<random>
#include<ctime>
using std::cout,std::vector;

// 元音/辅音的 两个 bool 参数之前设了默认值，由于后带了一个 singer，暂时没法用了，只能根据默认值重传一下。

void realizeDefaultADSR(int pos,int totsamp,const Note &note,vector<double> &buff,Singer &singer){
    for(int i=0;i<totsamp;i++){
        buff[pos+i]*=singer.getADSR(note.lyr,floor(100.0*i/totsamp));
    }
    // 注意先做这个 ADSR 再调响度，为了比较方便调教，把元音的默认 ADSR 调得保守一点，形状如一个梯形
}
void effectParameters(int pos,int totsamp,const Note &note,AudioFile<double>::AudioBuffer &buff,Singer &singer){
    // 注意因为各种音的 Delay 不同，所以这里把 totsamp 传进来。
    // 提前做一下响度压缩。（考虑把这个做成可选项）
    for(int channel=0;channel<numChannels;channel++){
        // dynCompress(pos,totsamp,buff[channel]);
        realizeDefaultADSR(pos,totsamp,note,buff[channel],singer);
        updateDyn(pos,totsamp,note,buff[channel]);
        updateEQ(pos,totsamp,note,buff[channel]);
    }
}

void GenCon(int &pos,const Note &note,AudioFile<double>::AudioBuffer &buff,bool KeepFormant,Singer &singer){    // Defult: Keep... = false.
    if(DEBUG)cout<<"Generating note: "<<note.lyr<<std::endl;
    int totsamp=beatToSamples(note.beat,BPM,sampleRateGlobal)+conDelaySamples;
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
        FFTfilter(pos,totsamp,buff[channel],buff[channel],sampleRateGlobal,req,1);
   
    effectParameters(pos,totsamp,note,buff,singer);
        
    pos+=totsamp-conDelaySamples;
}
void GenVow(int &pos,const Note &note,AudioFile<double>::AudioBuffer &buff,bool updPos,bool NoFormant,Singer &singer){  // Defult: true, false
    if(DEBUG)cout<<"Generating note: "<<note.lyr<<std::endl;
    int totsamp=beatToSamples(note.beat,BPM,sampleRateGlobal)+vowDelaySamples;  // 合成后置衔接过渡
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
                    b=MathCalc::Mod((w-now*2.0*M_PI)*(i-1)/sampleRateGlobal+b,2.0*M_PI);  // 我就说这里 b 怎么一直是 0，原来是上面没 floor
                    // Note：观察波形可以知道，后面这段完全不匀称，显然是相位有了问题。
                    w=now*2.0*M_PI;
                }
            }
            for(int channel=0;channel<numChannels;channel++){
                buff[channel][to]+=k*sin((static_cast<double>(i)/sampleRateGlobal)*w+b);
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
            FFTfilter(pos,totsamp,buff[channel],buff[channel],sampleRateGlobal,req,1);
    };
    double freq=keyToFreq(note.key);
    double now=freq;
    if(note.lyr=="SP")  goto finished;
    for(int i=0;i<20 && now<=18000;i++){ 
        genAccordingFreq(now,1.0/20,i+1);    // 记得平均。
        now+=freq;  
    }
    
    FormantEQ();
    effectParameters(pos,totsamp,note,buff,singer);

    finished:
    if(DEBUG)cout<<"Generated!"<<std::endl;
    if(updPos)
        pos+=totsamp-vowDelaySamples;
}
/*
上面的 FFT 使用，其实是一个误操作。当时写的时候 MathCalc 出了问题。所以改了。
上面事实是不需要 FFT 的。
一个事实是，如果能够把这些个东西统一处理就不会爆炸。得到音量比例，统一缩放就好了。
*/

/*

GenMix
old version:

auto genAccordingFreq=[&](double f0,double k){  // NOTE: 考虑修改这一部分。
    double k1=singer.getCoefficient(X.v1,f0),k2=singer.getCoefficient(X.v2,f0);
    double beg=X.tm1*(1-X.mixK),ed=X.tm1+X.tm2*X.mixK,nowK;
    beg*=totsamp,ed*=totsamp;
    for(int i=0;i<totsamp;i++){
        int to=pos+i;
        if(i<beg)   nowK=k1;
        else if(i<ed)   nowK=(k2-k1)/(ed-beg)*(i-beg)+k1;
        else    nowK=k2;
        for(int channel=0;channel<numChannels;channel++){
            buff[channel][to]+=k*nowK*sin((static_cast<double>(i)/sampleRateGlobal)*f0*2.0*M_PI);
        }
    }
};

double freq=keyToFreq(note.key);
double now=freq;
for(int i=0;i<20;i++){ 
    genAccordingFreq(now,1.0/20);    // 记得平均。
    now+=freq;
}

*/
void GenMix(int &pos,const Note &note,AudioFile<double>::AudioBuffer &buff,Singer &singer){
    if(DEBUG)cout<<"Generating note: "<<note.lyr<<std::endl;
    int totsamp=beatToSamples(note.beat,BPM,sampleRateGlobal)+vowDelaySamples;  // 合成后置衔接过渡
    auto &X=singer.Mix[note.lyr];
    
    auto &pit=note.pit;
    auto genframe_test=[&](){
        double f0=keyToFreq(note.key);
        double beg=X.tm1*(1-X.mixK),ed=X.tm1+X.tm2*X.mixK;
        beg*=totsamp,ed*=totsamp;   // added: 复制掉了一行。
        vector<double> w(20),b(20);
        for(int lay=0;lay<20;lay++)
            w[lay]=(pit.empty()?f0:pit[0].second+f0) * (lay+1)*2.0*M_PI,b[lay]=0,cout<<"wlay: "<<w[lay]<<std::endl;   // fixed: 乘号放外面记得。这里和上面 GenVow 不太一样了。
            // f0 是真的基频了。
        //
        // if(DEBUG)   cout<<"tot: "<<totsamp<<std::endl;
        for(int i=0,j=0,pitsz=pit.size();i<totsamp;i++){
            int to=pos+i;
            while(j+1<pitsz && i>=pit[j+1].first*totsamp)
                ++j;
            for(int lay=0;lay<20;lay++){
                double f_lay=f0*(lay+1),frq=f_lay;
                if(pitsz && i>0){        
                    frq=MathCalc::getLinearVal(pit[j].first*totsamp,pit[j].second,
                                            pit[j+1].first*totsamp,pit[j+1].second,i)*(lay+1)+f_lay;
                    if(fabs(w[lay]-frq*2.0*M_PI)>eps*40){   // 规避一些精度问题
                        b[lay]=MathCalc::Mod((w[lay]-frq*2.0*M_PI)*(i-1)/sampleRateGlobal+b[lay],2.0*M_PI);
                        w[lay]=frq*2.0*M_PI;
                    }
                }
                double k1=singer.getCoefficient(X.v1,f_lay),k2=singer.getCoefficient(X.v2,f_lay),nowK;   // 这里我也不知道为什么使用 f_lay 效果显著更好，本来应该是 frq 的啊。
                // 大概是因为共振峰的拟合太差了，这里直接上 f_lay 的固定比例还更好一点。
                // if(note.lyr=="aim")std::cerr<<k1<<" "<<k2<<std::endl;
                if(i<beg)   nowK=k1;
                else if(i<ed)   nowK=(k2-k1)/(ed-beg)*(i-beg)+k1;
                else    nowK=k2;
                nowK*=0.06;
                // if(DEBUG && i%10000==0)
                    // cout<<"nowK: "<<nowK<<std::endl;
                for(int channel=0;channel<numChannels;channel++){
                    buff[channel][to]+=nowK*sin((static_cast<double>(i)/sampleRateGlobal)*w[lay]+b[lay]);
                }
            }
            // NOTE：如果在这里强制使用之前的音高不变的合成方法就会爆炸，只能上 EQ 了。
        }
    };
    genframe_test();
    
    effectParameters(pos,totsamp,note,buff,singer);
    
    if(DEBUG)cout<<"Generated!"<<std::endl;
    pos+=totsamp-vowDelaySamples;
}

#endif