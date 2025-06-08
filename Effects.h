#ifndef EFFECT_H
#define EFFECT_H 1

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

#endif