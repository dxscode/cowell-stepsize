#ifndef Baseoper_h
#define Baseoper_h

class TVector
{ public:
     TVector();
     double& operator[](int i);
     TVector& operator=(TVector r);
     TVector& operator=(const double(&r)[3]);
     TVector operator+(TVector r);
     TVector operator-(TVector r);
     bool operator>(TVector r);
     double operator*(TVector r);
     TVector operator*(double s);
     TVector operator/(double s);
     TVector& VZero();
     TVector& I();
     double DistanceTo(TVector& r);
     double Length();
     TVector Cross(TVector r);
  private:
     double v[3];
};

class TMatrix
{ public:
    TMatrix();
    TMatrix(const TMatrix &r);
    TMatrix(double r[3][3]);
    double& operator()(int i,int j);
    TMatrix& operator=(TMatrix r);
    TMatrix operator+(TMatrix r);
    TMatrix operator-(TMatrix r);
    TVector operator*(TVector r);
    TMatrix operator*(double s);
    TMatrix MTransp();
    void I();
  private:
    double m[3][3];
};




TMatrix Rotate(double alpha,int axis);
double XyToRs(double x,double y,double &r,double &s);
double XyzToRls(TVector x,double &r,double &l,double &s);
double AdjustTo0_2PI(double s);
double AdjustTo_PItPI(double &s);
void ElemToCor(double mu,double dt,double elem[],TVector &r,TVector &v);
void CorToElem(double mu,double dt,TVector &r,TVector &v,double elem[]);
bool CorToCor(double CenterGm,double dt,TVector x0,TVector v0,TVector &x,TVector &v);

double uniform(double s,double e);
double Gaussian(double u,double sigma);
#endif
