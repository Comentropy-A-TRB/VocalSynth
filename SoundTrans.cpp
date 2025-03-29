#include<cstdio>
#include<cmath>
#include<cassert>
#include<string>
#include<vector>
#include<utility>
#include<tuple>
#include "./fftw/fftw3.h"

using std::string,std::vector;
/*
g++ .\SoundTrans.cpp main.cpp -o main -O2 -Wall -std=c++17 -L"./fexe: cannot find -libfftw3l-3ftw/" -llibfftw3-3 -llibfftw3f-3 -llibfftw3l-3
*/
const double A1=55.0;
const double pitD=pow(2,1.0/12);
const int tabDelta[12]={0,2,-9,-7,-5,-4,-2,0};

double keyToFreq(const string &s){    // A1 -> 55Hz
    int up=0;
    if(s[0]=='#')   up=1;
    assert((int)s.size()==2+up && isalpha(s[up]) && isdigit(s[up+1]));
    return A1*pow(pitD,tabDelta[s[up]-'A']+up)*pow(2,s[up+1]-'1');
}
double beatToTime(double beat,int BPM){
    return beat*60/BPM;
}
double timeToSamples(double t,double sampleRate){
    return t*sampleRate;
}
double beatToSamples(double beat,int BPM,double sampleRate){
    return timeToSamples(beatToTime(beat,BPM),sampleRate);
}
/**
 * @Description : 使用FFT进行滤波
 *                使用示例：
 *                  原始采样频率为100kHz，采集了10000个点，保存为单精度浮点数。滤除其中20kHz~30kHz的频率
 *                  fft_filter_f(10000, in, out, 100e3, 20e3, 30e3)
 * @Input       : n             输入数据个数
 *                in            输入数据
 *                              in和out指向不同位置，不改变输入
 *                              in和out指向相同位置，输出将覆盖输入
 *                sample_rate   原始数据采样频率 Hz
 *                freq_start    需过滤的起始频率 Hz
 *                freq_end      需过滤的截止频率 Hz
 *                              若freq_end大于采样频率的50%，则将滤除大于freq_start的所有高频信号
 * @Output      : out           输出数据
 *                              in和out指向不同位置，不改变输入
 *                              in和out指向相同位置，输出将覆盖输入
 * @Return      : 无
*/

void FFTfilter(int startPos,int n, const vector<double> &in, vector<double> &out,double sampleRate,vector<std::tuple<double,double,double> > &req,int opt=0){
    // opt==0 表示赋值， opt==1 表示在基础上 EQ
    int i, begin, end;
    double *real;
    fftw_complex *complex;
    fftw_plan plan;

    // fftw的内存分配方式和malloc类似，但使用SIMD（单指令多数据流）时，fftw_alloc会将数组以更高效的方式对齐
    real = fftw_alloc_real(n);
    complex = fftw_alloc_complex(n/2+1);    // 实际只会用到(n/2)+1个complex对象

    // Step1：FFT实现时域到频域的转换
    plan = fftw_plan_dft_r2c_1d(n, real, complex, FFTW_ESTIMATE);
    for (i = 0; i < n; i++)
    {
        real[i] = in[i + startPos];
    }

    // 对长度为n的实数进行FFT，输出的长度为(n/2)-1的复数
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Step2：计算需滤波的频率在频域数组中的下标
    for(auto &[freqStart,freqEnd,loud]:req){
        begin = (int)((freqStart / sampleRate) * n);
        end = (int)((freqEnd / sampleRate) * n);
        end = end < (n / 2 + 1) ? end : (n / 2 + 1);
        for (i = begin; i < end; i++)
        {
            // 对应的频率分量置为0，即去除该频率
            // upd: 改为乘以一个系数
            if(opt==0){
                complex[i][0] = loud;
                complex[i][1] = loud;
            }else{
                complex[i][0] *= loud;
                complex[i][1] *= loud;
            }
        }
    }
    // Step3：IFFT实现频域到时域的转换
    // 使用FFTW_ESTIMATE构建plan不会破坏输入数据
    plan = fftw_plan_dft_c2r_1d(n, complex, real, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Step4：计算滤波后的时域值
    for (i = 0; i < n; i++)
    {
        // 需除以数据个数，得到滤波后的实数
        out[i+startPos] = real[i] / n;
    }

    fftw_free(real);
    fftw_free(complex);
}