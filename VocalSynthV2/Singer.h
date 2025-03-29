#ifndef SINGER_H
#define SINGER_H 1

#include<vector>
#include<map>
#include<string>
#include<cassert>
#include<cmath>
#include<algorithm>
#include<utility>
using std::map,std::vector,std::string,std::pair;
const double eps=1e-6;
using constr=const string;
struct Singer{
    // Singer 内主要记默认参数，尤其共振峰比较重要。
    struct Formant{
        // 2 -> 直接按泛音生成一点，不用共振峰了。
        vector<int> F;
        vector<double> loud;
    };  // upd: 若干共振峰，数字在文件内部给出

    struct ADSR{    // begin from 0 ! ! !
        vector<int> T,loud;    // precent ;
        ADSR():T(4),loud(4){}
    };

    struct MixVow{  // 对于交变参数单独处理。
        string v1,v2;
        double tm1,tm2,mixK;
    };
    // 0 -> 辅音; 1 -> vowels, 2 -> 浊辅音，3-> 交变双元音
    map<string,int> typ;
    map<string,MixVow> Mix;
    map<string,Formant> Form;
    map<string,ADSR> timeLoud;
    double getQuadraticCurveVal(vector<double> X,vector<double> Y,double x){    // 拉格朗日插值三点式
        double ret=0;
        for(int i=0;i<3;i++){
            double now=Y[i];
            for(int j=0;j<3;j++)
                if(i!=j)    now=now*(x-X[j])/(X[i]-X[j]);
            ret+=now;
        }
        return ret;
    }
    double getQuadraticCurveVal(double midx,double midy,double nx,double ny,double x){  // 通过顶点式
        double K=(ny-midy)/pow(nx-midx,2);
        return K*pow(x-midx,2)+midy;
    }
    double getMediumX(double l,double r){
        if(r-l<=500)    return (l+r)/2;
        if(r-l<=1000)   return (l+2*r)/3;
        return (l+3*r)/4;
    }
    double getMediumVal(double l,double r){ // 输入两个共振峰返回某个中间点的二次函数 y 坐标
        if(r-l<=500) return 0.5;
        if(r-l<=1000)   return 0.2;
        return 0.1;
    }
    double getCoefficient(constr &lyr,double freq){
        auto &fx=Form[lyr];
        if(freq<=fx.F[1])
            return getQuadraticCurveVal(0.0,0.1,fx.F[1],fx.loud[1],freq);
        for(int i=1;i<=2;i++)
            if(freq<=fx.F[i+1])
                return getQuadraticCurveVal({(double)fx.F[i],(double)fx.F[i+1],getMediumX(fx.F[i],fx.F[i+1])},
                                        {fx.loud[i],fx.loud[i+1],getMediumVal(fx.F[i],fx.F[i+1])},freq);
        return 0.03;
    }
    double getConCoeffic(constr &lyr,double freq){
        auto &fx=Form[lyr];  // 辅音的第二共振峰填入一个类周期频率。
        auto dist=[&](int num,double k){
            return fx.loud[num]-k*fabs(fx.F[num]-freq);
        };
        if(freq<=fx.F[1])   
            return std::max(eps,dist(1,0.0004));
        if(freq>=fx.F[3])
            return std::max(eps,dist(3,0.0004));
        return std::max({eps,dist(1,0.0005),dist(2,0.0005),dist(3,0.0005)});
    }
    double getADSR(constr &lyr,int perc){
        auto &fx=timeLoud[lyr];
        auto getPointADSRVal=[&](int sx,int sy,int ex,int ey){
            if(sx-ex==0)    return (double)sy;
            return 1.0*(sy-ey)/(sx-ex)*(perc-sx)+sy;
        };
        for(int i=0,sz=fx.T.size();i+1<sz;i++){
            if(fx.T[i]<=perc && perc<=fx.T[i+1])    // 记得乘一个 1/100
                return 0.01*getPointADSRVal(fx.T[i],fx.loud[i],fx.T[i+1],fx.loud[i+1]);
        }
        assert(false);
        return 0;
    }
    void setType(constr &lyr,int _typ){
        typ[lyr]=_typ;
    }
    void setFormant(constr &lyr,vector<int> &_F,vector<double> &_loud){
        Formant &X=Form[lyr];
        X.F=_F,X.loud=_loud;
    }
    void setADSR(constr &lyr,vector<int> &_T,vector<int> &_loud){
        ADSR &X=timeLoud[lyr];
        X.T=_T,X.loud=_loud;
    }
    void setMixVow(constr &Vow,constr &v1,constr &v2,vector<double> &K){
        MixVow &X=Mix[Vow];
        X={v1,v2,K[0],K[1],K[2]};
    }
    int getType(constr &lyr){
        return typ[lyr];
    }
};
#endif 
