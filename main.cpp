
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TBody.h"

using namespace std;


const double AU=1.495978706910000E+08;
const double SunGm=2.959122082855911025e-04;


double _mu=2.959122082855911025e-04;
double _InitialElem[6]={1.5,0.01,1,1,1,1};
double _TC=0.00000000000000000000000000;

double _EndT=3652500;

double _DoubleCnt, _HalfCnt;

TBody _TBD;

TVector CalcuAcc(TVector R)
{ TVector A;
  double rr=R.Length();

  A=R*(-_mu/(rr*rr*rr));

  return A;
}


double _FT_back[3][15];

void CopyFT(int endk)
{ for(int j=0;j<3;j++) for(int k=-endk;k<=endk;k++) _FT_back[j][k+7]=_TBD.f[j][k+7];
}

double CmpFT(int endk)
{ double sum=0.00000000000000000000000000000;
  for(int j=0;j<3;j++) for(int k=-endk;k<=endk;k++)
    { if(_FT_back[j][k+7]==0.00000000000000000000000)
         sum+=fabs(_TBD.f[j][k+7]);
      else
         sum+=fabs((_TBD.f[j][k+7]-_FT_back[j][k+7])/_FT_back[j][k+7]);
    }
  return sum;
}

void Initial(int order,double h)
{ TVector R0,V0,A,R,V,A0,R1,V1,A1,Rf1,Vf1,Af1,Af1back;
  double h2,epslong;
  int ord;

  _TBD.h=h;         _DoubleCnt=0; _HalfCnt=0;
  h2=h;
  epslong=1.0e-25;
  _TC=0.0000000000000000000000;

  ElemToCor(_mu,0.0,_InitialElem,R0,V0);
  A0=CalcuAcc(R0);
  A1=A0; Af1=A0;

  _TBD.WriteF(-1,Af1);   _TBD.WriteF(0,A0);    _TBD.WriteF(1,A1);


  for(int ord=2;ord<=12;ord+=2)
    { do{ CopyFT(ord/2);

          _TBD.CalcuCentralSf3(ord,R0,V0);
          for(int k=1;k<=ord/2;k++)
            { _TBD.CalcuIntXv(k,ord,k,R,V);   A=CalcuAcc(R);
              _TBD.WriteF(k,A);
              _TBD.CalcuCentralSf3(ord,R0,V0);
              _TBD.CalcuSfAccCentral3(ord/2+1);
            }

          for(int k=-1;k>=-ord/2;k--)
            { _TBD.CalcuIntXv(k,ord,k,R,V);   A=CalcuAcc(R);
              _TBD.WriteF(k,A);
              _TBD.CalcuCentralSf3(ord,R0,V0);
              _TBD.CalcuSfAccCentral3(ord/2+1);
            }
        }while(CmpFT(ord/2)>epslong);

      _TBD.CalcuCentralSf3(ord,R0,V0);
      _TBD.CalcuSfAccCentral3(ord/2+1);

      _TBD.CalcuIntXv(ord/2+1,ord,ord/2+1,R,V);   A=CalcuAcc(R);
      _TBD.WriteF(ord/2+1,A);
      _TBD.CalcuIntXv(ord/2+1,ord,ord/2,R,V);   A=CalcuAcc(R);
      _TBD.WriteF(ord/2+1,A);

      _TBD.CalcuIntXv(-ord/2-1,ord,-ord/2-1,R,V);   A=CalcuAcc(R);
      _TBD.WriteF(-ord/2-1,A);
      _TBD.CalcuIntXv(-ord/2-1,ord,-ord/2,R,V);   A=CalcuAcc(R);
      _TBD.WriteF(-ord/2-1,A);
    }

  _TBD.CalcuCentralSf3(order,R0,V0);
  _TBD.CalcuSfAccCentral3(7);
}

void sbs_ord(int ord)
{ TVector R,V,A;

  _TBD.sbs(1,7);
  _TBD.CalcuIntXv(7,ord,ord/2+1,R,V);
  A=CalcuAcc(R);
  _TBD.WriteF(7,A);

  _TBD.MoveForwardOnce();
  _TC+=_TBD.h;

  _TBD.CalcuIntXv(6,ord,ord/2,R,V);
  A=CalcuAcc(R);
  _TBD.WriteF(6,A);

}

void HalfPace()
{ TVector R,V,A;

  for(int k=-5;k<6;k+=2)
    { long double th=k;
      _TBD.CalcuFloatXv(th/2,R,V);
      A=CalcuAcc(R);
      for(int m=0;m<3;m++) _TBD.fb[m][k+7]=(A[m]*_TBD.h*_TBD.h)/4;
    }
  for(int k=-6;k<7;k+=2)
      for(int m=0;m<3;m++) _TBD.fb[m][k+7]=_TBD.f[m][k/2+7]/4;

  _TBD.CalcuIntXv(0,R,V);
  for(int i=1;i<14;i++)  for(int m=0;m<3;m++) _TBD.f[m][i]=_TBD.fb[m][i];
  _TBD.h/=2;
  _TBD.CalcuCentralSf3(12,R,V);
  _TBD.CalcuSfAccCentral3(7);

  _HalfCnt+=1;
}

void DoublePace(long double tc,TVector Rc,TVector Vc)
{ for(int i=1;i<14;i++)  for(int m=0;m<3;m++) _TBD.f[m][i]=_TBD.fb[m][i];
  _TBD.h*=2;
  _TBD.CalcuCentralSf3(12,Rc,Vc);
  _TBD.CalcuSfAccCentral3(7);
  _TC=tc;

  _DoubleCnt+=1;
}




int main()
{ _mu=SunGm;
  _InitialElem[0]=1;  _InitialElem[1]=0.01; for(int i=2;i<6;i++) _InitialElem[i]=1;
  //-*+*+*+*+*+**+*+*-*+*+*+*+*+**+*+*-*+*+*+*+*+**+*+*-*+*+*+*+*+**+*+*-*+*+*+

  ofstream out("C24.txt");
  int ord=12,otc=0,pc;
  TVector R,V,Ra,Va,Rc,Vc;
  double maxe,d,dv,dt,tcback;

  dt=_EndT/10000;

  Initial(ord,4);
  maxe=-1;  otc=0;   pc=0;
  while(_TC<_EndT)
    { if(_TC+0.1>otc)
        { _TBD.CalcuFloatXv((otc-_TC)/_TBD.h,R,V);
          ElemToCor(_mu,otc,_InitialElem,Ra,Va);
          d=R.DistanceTo(Ra);  dv=V.DistanceTo(Va);
          if(d>maxe) maxe=d;
          out<<otc/365.25<<"  "<<d<<"  "<<dv<<endl;

          otc+=dt;
        }

      if(pc%15==0) HalfPace();
      if(pc%15==2) for(int m=0;m<3;m++) for(int k=-6;k<=6;k+=2) _TBD.fb[m][(k-6)/2+7]=_TBD.f[m][k+7]*4;
      if(pc%15==4) for(int m=0;m<3;m++) _TBD.fb[m][1+7]=_TBD.f[m][6+7]*4;
      if(pc%15==6) for(int m=0;m<3;m++) _TBD.fb[m][2+7]=_TBD.f[m][6+7]*4;
      if(pc%15==8) { for(int m=0;m<3;m++) _TBD.fb[m][3+7]=_TBD.f[m][6+7]*4;
                      _TBD.CalcuIntXv(0,Rc,Vc);  tcback=_TC;
                   }
      if(pc%15==10) for(int m=0;m<3;m++) _TBD.fb[m][4+7]=_TBD.f[m][6+7]*4;
      if(pc%15==12) for(int m=0;m<3;m++) _TBD.fb[m][5+7]=_TBD.f[m][6+7]*4;
      if(pc%15==14) { for(int m=0;m<3;m++) _TBD.fb[m][6+7]=_TBD.f[m][6+7]*4;
                      DoublePace(tcback,Rc,Vc);
                    }

      _TBD.CalcuIntXv(0,R,V);
      ElemToCor(_mu,_TC,_InitialElem,Ra,Va);
      d=R.DistanceTo(Ra);
      if(d>maxe) maxe=d;

      sbs_ord(ord);    pc++;
    }
  cout<<_DoubleCnt<<", "<<_HalfCnt<<endl;
  
  return 0;
}
//---------------------------------------------------------------------------




