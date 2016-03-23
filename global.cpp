#include "global.h"
//---------------------------------------------------------------------------
string global::getFileName(string pathstring) {
  boost::filesystem::path p(pathstring);
  return p.filename().c_str();
  }
//---------------------------------------------------------------------------
void CALC::sran1(long seedvalue) {
  // Initialize with negative number
  if (rseed<0)
    rseed=seedvalue;
  }
//---------------------------------------------------------------------------
double CALC::ran1() {
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (rseed<=0 || !iy) {
    if (-rseed<1)
      rseed=1;
    else
      rseed=-rseed;
    for (j=NTAB+7; j>=0; j--) {
      k=rseed/IQ;
      rseed=IA*(rseed-k*IQ)-IR*k;
      if (rseed<0)
        rseed+=IM;
      if (j<NTAB)
        iv[j]=rseed;
      }
    iy=iv[0];
    }
  k=(rseed)/IQ;
  rseed=IA*(rseed-k*IQ)-IR*k;
  if (rseed<0)
    rseed+=IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=rseed;
  temp=AM*iy;
  if (temp>RNMX)
    return RNMX;
  return temp;
  }
//---------------------------------------------------------------------------
double CALC::Chi2(double x, int n) {
  double p,t;
  int a,k;

  if (n==0)
    return 1;
  p=exp(-0.5*x);
  if ((n%2)==1)
    p*=sqrt(2*x/M_PI);
  for (k=n; k>=2; k-=2)
    p*=x/(double)k;
  t=p;
  a=n;
  while (t>0.000001*p) {
    a+=2;
    t*=x/(double)a;
    p+=t;
    }
  return 1-(p<0?0:p>1?1:p);
  }
//---------------------------------------------------------------------------
