#ifndef SOUNDTRANS_H
#define SOUNDTRANS_H 1

#include<string>
#include<tuple>
#include<vector>
double keyToFreq(const std::string &s);
double beatToTime(double beat,int BPM);
double timeToSamples(double t,double sampleRate);
double beatToSamples(double beat,int BPM,double sampleRate);
void FFTfilter(int startPos,int n, const std::vector<double> &in, std::vector<double> &out,double sampleRate,std::vector<std::tuple<double,double,double> > &req,int opt=0);

#endif 
