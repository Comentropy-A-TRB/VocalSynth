#ifndef MATH_CALC_H
#define MATH_CALC_H

#include "Setting.h"
#include<cmath>
#include<vector>
#include<numeric>

namespace MathCalc{
    using std::vector;
    double getLinearVal(double x0,double y0,double x1,double y1,double x){
        if(fabs(x0-x1)<eps) return (y0+y1)/2;
        return 1.0*(y0-y1)/(x0-x1)*(x-x0)+y0;
    }
    double Mod(double x,double y){
        return x-floor(x/y)*y;  // 笑死我了没写 floor.
    }
    void vectorScaling(vector<double> &vec,double unit=1){
        double ratio=std::accumulate(vec.begin(),vec.end(),0)/unit;
        for(auto &i:vec)    i/=ratio;
    }
}

#endif