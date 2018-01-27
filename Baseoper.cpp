#include "Baseoper.h"
#include <math.h>
#include <stdlib.h> 

using namespace std;

TVector::TVector()
        { for(int i=0;i<3;i++)  v[i]=0.0;
        }
     double& TVector::operator[](int i)
        { return v[i];
        }
     TVector& TVector::operator=(TVector r)
        { for(int i=0;i<3;i++)
            v[i]=r[i];
          return *this;
        }
     TVector& TVector::operator=(const double(&r)[3])
        { for(int i=0;i<3;i++)
            v[i]=r[i];
           return *this;
        }
     TVector TVector::operator+(TVector r)
        { TVector p;
          for(int i=0;i<3;i++)  p[i]=v[i]+r[i];
          return p;
        }
     TVector TVector::operator-(TVector r)
        { TVector p;
          for(int i=0;i<3;i++)  p[i]=v[i]-r[i];
          return p;
        }
     double TVector::operator*(TVector r)
        { double s=0.0;
          for(int i=0;i<3;i++)  s+=v[i]*r[i];
          return s;
        }
     TVector TVector::operator*(double s)
        { TVector d;
          for(int i=0;i<3;i++)  d[i]=v[i]*s;
          return d;
        }
     bool TVector::operator>(TVector r)
        { double r1,r2;
          r1=(*this)*(*this);
          r2=r*r;
          if(r1>r2) return true;
          return false;
        }
     TVector& TVector::VZero()
        { for(int i=0;i<3;i++)  v[i]=0.0;
          return *this;
        }
     double TVector::DistanceTo(TVector& r)
        { double dis=sqrt((v[0]-r[0])*(v[0]-r[0])+(v[1]-r[1])*(v[1]-r[1])+(v[2]-r[2])*(v[2]-r[2]));
          return dis;
        }
     double TVector::Length()
        { return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        }
     TVector& TVector::I()
        { double r=this->Length();
          for(int i=0;i<3;i++)  v[i]/=r;
          return *this;
        }
     TVector TVector::Cross(TVector r)
       { TVector p;
         p[0]=v[1]*r[2]-v[2]*r[1];
         p[1]=v[2]*r[0]-v[0]*r[2];
         p[2]=v[0]*r[1]-v[1]*r[0];
         return p;
       }
    TVector TVector::operator/(double s)
        { TVector d;
          for(int i=0;i<3;i++)  d[i]=v[i]/s;
          return d;
        }

    TMatrix::TMatrix()
        { for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
               m[i][j]=0.0;
        }
    TMatrix::TMatrix(const TMatrix &r)
        { for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
               m[i][j]=r.m[i][j];
        }
    TMatrix::TMatrix(double r[3][3])
        { for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
               m[i][j]=r[i][j];
        }
    double& TMatrix::operator()(int i,int j)
        { return m[i][j];
        }
    TMatrix& TMatrix::operator=(TMatrix r)
        { for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
               m[i][j]=r(i,j);
          return *this;
        }
    TMatrix TMatrix::operator+(TMatrix r)
        { TMatrix p;
          for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
               p(i,j)=m[i][j]+r(i,j);
          return p;
        }
    TMatrix TMatrix::operator-(TMatrix r)
        { TMatrix p;
          for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
               p(i,j)=m[i][j]-r(i,j);
          return p;
        }
    TVector TMatrix::operator*(TVector r)
        { TVector p;
          for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
               p[i]+=m[i][j]*r[j];
          return p;
        }
    TMatrix TMatrix::operator*(double s)
        { TMatrix p;
          for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
               p(i,j)=m[i][j]*s;
          return p;
        }
    TMatrix TMatrix::MTransp()
        { TMatrix p;
          for(int i=0;i<3;i++)
             for(int j=0;j<3;j++)
                p(i,j)=m[j][i];
          return p;
        }
    void TMatrix::I()
        { for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
               m[i][j]=(i==j)?1.0L:0.0L;
        }





TMatrix Rotate(double alpha,int axis)
{  double xsin=sin(alpha);
   double xcos=cos(alpha);
   int i=axis,j=axis;
   TMatrix r;
   r(i-1,j-1)=1.0;
   if(++i>3) i=i-3;
   if(++j>3) j=j-3;
   r(i-1,j-1)=xcos;
   if(++j>3) j=j-3;
   r(i-1,j-1)=xsin;
   r(j-1,i-1)=-xsin;
   if(++i>3) i=i-3;
   r(i-1,j-1)=xcos;
   return r;
}
double XyToRs(double x,double y,double &r,double &s)
{ r=sqrt(x*x+y*y);
  s=y>0?acos(x/r):2*M_PI-acos(x/r);
  return r;
}
double XyzToRls(TVector x,double &r,double &l,double &s)
{ double r1=XyToRs(x[0],x[1],r,l);
  XyToRs(r1,x[2],r,s);
  s=s>M_PI?s-2*M_PI:s;
  return r;
}
double AdjustTo0_2PI(double s)
{ double a,b,c,d;
  b=s/2/M_PI;
  a=modf(b,&c);
  d=s-2*M_PI*c;
  d=d<0?d+2*M_PI:d;
  return d;
}
double AdjustTo_PItPI(double &s)
{ while(s<M_PI)
   { s+=2*M_PI;
   }
  while(s>M_PI)
   { s-=2*M_PI;
   }
  return s;
}

double SolveKepler(double e,double Mf)
{ double M,E,c,max,min,mid,fmin,fmid,fmax;
  int cnt;
  bool ff;

  M=AdjustTo0_2PI(Mf);
  if(M==0||M==M_PI||M==M_PI*2) return M;

  ff=M>3.1415926535897932384626433832795L;
  if(ff)
    M=6.283185307179586476925286766559L-M;

  E=M;    cnt=0;
  do{ c=(M-E+e*sin(E))/(1-e*cos(E));
      E+=c;
      cnt++;
    }while(fabs(c)>1.0e-15&&cnt<20);

  if(cnt<20)
     return ff?-E:E;

  max=M_PI>(M+e)?M+e:M_PI;
  min=M;

  fmax=max-e*sin(max)-M;     // if(fabs(fmax)<1e-15) return max;
  fmin=min-e*sin(min)-M;     // if(fabs(fmin)<1e-15) return min;

  do{ mid=(min+max)/2;
      fmid=mid-e*sin(mid)-M;
      if(fabs(fmid)<1e-15)
         return ff?-mid:mid;

      if(fmid*fmax>0) {max=mid; fmax=fmid;}
      else            {min=mid; fmin=fmid;}
    }while(max-min>1.0e-15);

  return ff?-mid:mid;
}

void ElemToCor(double mu,double dt,double elem[],TVector &r,TVector &v)
{ TVector P,Q;
  P[0]=1.0; Q[1]=1.0;
  double dw=elem[3],i=elem[4],w=elem[2];
  P=Rotate(-dw,3)*(Rotate(-i,1)*(Rotate(-w,3)*P));            //elem[0]  a
  Q=Rotate(-dw,3)*(Rotate(-i,1)*(Rotate(-w,3)*Q));            //elem[1]  e
  double n,M,E,p1,p2,c;                                  //elem[2]  omigar
  n=sqrt(mu/(elem[0]*elem[0]*elem[0]));                      //elem[3]  Node
  M=n*dt+elem[5];                                             //elem[4]  i
  M=AdjustTo0_2PI(M);                                         //elem[5]  M0

  E=SolveKepler(elem[1],M);
                           /*               
  E=M;
  do{ c=(M-E+elem[1]*sin(E))/(1-elem[1]*cos(E));
      E+=c;
    }while(fabs(c)>1.0e-14); /* */

  p1=elem[0]*(cos(E)-elem[1]);
  p2=elem[0]*sqrt(1-elem[1]*elem[1])*sin(E);
  r=P*p1+Q*p2;
  p1=-sin(E)*elem[0]*elem[0]*n/sqrt(r*r);
  p2=elem[0]*elem[0]*n*sqrt(1-elem[1]*elem[1])*cos(E)/sqrt(r*r);
  v=P*p1+Q*p2;
 // r=Rotate(-23.43928108*M_PI/180.0,1)*r;
 // v=Rotate(-23.43928108*M_PI/180.0,1)*v;

}

void CorToElem(double mu,double dt,TVector &r,TVector &v,double elem[])
{ double SumGm,rr,vr2,a,n,eccen,m0,node,omigar,i,f,
              ecosE,esinE,E,k1,k2,temp,sini;
  TVector R,V,P,Q;
  R=r;
  V=v;
  SumGm=mu;
  rr=sqrt(R*R);
  vr2=V*V;
  a=rr/(2.0-rr*vr2/SumGm);
  n=sqrt(SumGm/(a*a*a));
  ecosE=1.0-rr/a;
  esinE=(R*V)/(a*a*n);
  XyToRs(ecosE,esinE,eccen,E);
  m0=E-esinE-n*dt;
  k1=(ecosE-eccen*eccen)/(1-ecosE)/eccen;                 //-------
  k2=sqrt(1-eccen*eccen)*esinE/(1-ecosE)/eccen;          // f
  XyToRs(k1,k2,temp,f);                       // true anomaly
  k1=cos(E)/rr;
  k2=-sin(E)/(a*n);
  P=R*k1+V*k2;
  k1=sin(E)/(rr*sqrt(1-eccen*eccen));
  k2=(cos(E)-eccen)/(a*n*sqrt(1-eccen*eccen));
  Q=R*k1+V*k2;
  //---------------------------------
  XyToRs(Q[2],P[2],sini,omigar);                   // omigar
  k1=(P[1]*sin(omigar)+Q[1]*cos(omigar));
  k2=-(P[0]*sin(omigar)+Q[0]*cos(omigar));
  XyToRs(k1,k2,temp,node);                      // node
  XyToRs(k1/cos(node),sini,temp,i);
  //---------------------------------
  elem[0]=a;
  elem[1]=eccen;
  elem[2]=omigar;
  elem[3]=node;
  elem[4]=i;
  elem[5]=AdjustTo0_2PI(m0);
}

bool CorToCor(double CenterGm,double dt,TVector x0,TVector v0,TVector &x,TVector &v)
        { double SumGm,rr0,v2,a,n,dM,c,ecosE0,esinE0,eccen,dE,A,B,C,D,rr;
          SumGm=CenterGm;
          rr0=sqrt(x0*x0);
          v2=v0*v0;
          a=rr0*SumGm/(2*SumGm-rr0*v2);
          int ddd=0;
          if(a>0)
             { n=sqrt(SumGm/(a*a*a));
               ecosE0=(a-rr0)/a;
               esinE0=(x0*v0)/(a*a*n);

               dM=AdjustTo0_2PI(dt*n);
               //-------------------------------------------------------------------
               if(dM>M_PI) dM-=2*M_PI;
               dM=dt*n;
               dE=dM;
               ddd=0;
               do{ ddd++;
                  c=(dM-dE+ecosE0*sin(dE)-esinE0*(1-cos(dE)))/(1-ecosE0*cos(dE)+esinE0*sin(dE));
                   dE+=c;
                 }while(fabs(c)>1.0e-10&&ddd<100000);
               if(ddd>=100000) return false;
               //-------------------------------------------------------------------
               A=1-a*(1-cos(dE))/rr0;
               B=(dM-(dE-sin(dE)))/n;
               x=x0*A+v0*B;
               rr=sqrt(x*x);
               C=-a*a*n*sin(dE)/(rr0*rr);
               D=1-a*(1-cos(dE))/rr;
               v=x0*C+v0*D;
             }
          else
             { a=-a;
               n=sqrt(SumGm/(a*a*a));
               ecosE0=(a-rr0)/a;
               esinE0=(x0*v0)/(a*a*n);
               double ee,temp;
               XyToRs(ecosE0,esinE0,ee,temp);
               dM=dt*n;
               //-------------------------------------------------------------------
              // if(dM>M_PI) dM-=2*M_PI;    
               dE=dM;
               c=(dM+dE-ee*sinh(dE))/(1.0+ee*cosh(dE));
               ddd=0;
               while(fabs(c)>1.0e-10&&ddd<100000)
                 { ddd++;
                   dE+=c;
                  c=(dM+dE-ee*sinh(dE))/(1.0+ee*cosh(dE));
                 }
               if(ddd>=100000) return false;
               A=1-a*(1-cosh(dE))/rr0;
               B=(dM-(dE-sinh(dE)))/n;
               x=x0*A+v0*B;
               rr=sqrt(x*x);
               C=-a*a*n*sinh(dE)/(rr0*rr);
               D=1+a*(1-cosh(dE))/rr;
               v=x0*C+v0*D;
             }
          return true;
        }

int CharClass(char c)
// 0: NoView
// 1: 0~9
// 2: .
// 3: +-
// 4: a-z
// 5: A-Z
// 6: (+Shift)
// 7: (No Shift)
{ if(c=='.') return 2;
  if(c>='0'&&c<='9') return 1;
  if(c=='+'||c=='-') return 3;
  if('a'<=c&&c<='z') return 4;
  if('A'<=c&&c<='Z') return 5;
  if(c=='~'||c=='!'||c=='@'||c=='#'||c=='$'||c=='%'||c=='^'||c=='&'||c=='*'||c=='('||c==')'||c=='_'||c=='{'||c=='}'||c=='|'||c==':'||c=='"'||c=='<'||c=='>'||c=='?')
    return 6;
  if(c=='`'||c=='['||c==']'||c=='\\'||c==';'||c=='\''||c==','||c=='/'||c==' '||c=='=')
    return 7;
  return 0;
}


double uniform(double s,double e)
{ double ran;
  int dd=0x7FFFU;
  ran=(double)(rand());
  ran/=(double)dd;
  ran*=e-s;
  return ran+s;
}


double Gaussian(double u,double sigma)
{ double ray[20];
  for(int i=0;i<20;i++) ray[i]=uniform(-sigma/2.58,sigma/2.58);
  double s=0.0;
  for(int i=0;i<20;i++) s+=ray[i];
  return s+u;
}
