#ifndef TBody_h
#define TBody_h
#include "Baseoper.h"

class TBody
{ public:
    double h;
    double f[3][15];                 // k  -7  ----> +7
    double sf[3][15][3];                // k  -7  ----> +7
    TBody();
    void CalcuFloatXv(double dtVh,TVector& X,TVector& V);  // 计算非整数坐标速度    -7.0 <= dtVh <= +7.0
    void CalcuIntXv(int k,TVector& X,TVector& V);              // 计算整数步点坐标速度     k  -7  ----> +7
    void CalcuIntXv(int k_12,int order,int k_ord,TVector& X,TVector& V);   // 计算整数步点坐标速度     k_12:[-7,+7];  order:所用阶数;  k_ord:[-ord/2,ord/2]
    void CalcuCentralSf3(int order,TVector X,TVector V);     // 根据X、V计算中间三个步点的和分值，order是阶数
    void CalcuSfAccCentral3(int endk);
    void CalcuSfAccR(int k_12,int ord,int k_ord,TVector X);
    void sbs(int dir,int k_next);
    void WriteF(int k,TVector A);
    void MoveForwardOnce();
    void MoveBackOnce();
    TVector XError_12();
};

#endif



