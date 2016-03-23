#include "logreg.h"
#include "global.h"
//---------------------------------------------------------------------------
LogReg::LogReg() {
  theta=tmptheta=z=sumofsquares=NULL;
  dimx=dimy=0;
  }
//---------------------------------------------------------------------------
void LogReg::createArrays(double *y,double **xo,int dimy,int dimx) {
  this->dimy=dimy;
  this->dimx=dimx+1;
  this->y=y;
  this->xo=xo;
  x=global::make2DArray(x,this->dimy,this->dimx);
  theta=new double[this->dimx];
  tmptheta=new double[this->dimx];
  z=new double[this->dimx];
  sumofsquares=new double[this->dimx];
  }
//---------------------------------------------------------------------------
void LogReg::clearArrays() {
  fill_n(theta,dimx,0);
  for (int y1=0; y1<dimy; y1++)
    x[y1][0]=1;
  }
//---------------------------------------------------------------------------
double LogReg::matrixMultiply(double *n1,double *n2) {
  double res;

  res=0;
  for (int i1=0; i1<dimx; i1++)
    res+=n1[i1]*n2[i1];
  return res;
  }
//---------------------------------------------------------------------------
double LogReg::sigmoid(double z) {
  return 1/(1+exp(-z));
  }
//---------------------------------------------------------------------------
void LogReg::normalize() {
  double *minx,*maxx;

  minx=new double[dimx-1];
  maxx=new double[dimx-1];
  for (int y1=0; y1<dimy; y1++)
    for (int x1=0; x1<(dimx-1); x1++) {
      if (y1==0)
        minx[x1]=maxx[x1]=xo[y1][x1];
      else {
        minx[x1]=min(minx[x1],xo[y1][x1]);
        maxx[x1]=max(minx[x1],xo[y1][x1]);
        }
      }
  for (int y1=0; y1<dimy; y1++)
    for (int x1=0; x1<(dimx-1); x1++)
      x[y1][x1+1]=(xo[y1][x1]-minx[x1])/(maxx[x1]-minx[x1]);
  delete minx;
  delete maxx;
  }
//---------------------------------------------------------------------------
double LogReg::calculateCost() {
  double J,h0;

  J=0;
  for (int y1=0; y1<dimy; y1++) {
    h0=sigmoid(matrixMultiply(x[y1],theta));
    J=J-y[y1]*log(h0)-(1-y[y1])*log(1-h0);
    }
  return J/dimy;
  }
//---------------------------------------------------------------------------
void LogReg::calculateTheta() {
  double h0;

  fill_n(tmptheta,dimx,0);
  for (int y1=0; y1<dimy; y1++) {
    h0=sigmoid(matrixMultiply(x[y1],theta))-y[y1];
    for (int x1=0; x1<dimx; x1++) {
      tmptheta[x1]=tmptheta[x1]+h0*x[y1][x1];
      }
    }
  for (int x1=0; x1<dimx; x1++)
    tmptheta[x1]/=dimy;
  }
//---------------------------------------------------------------------------
void LogReg::calculateZ() {
  for (int x1=0; x1<dimx; x1++)
    for (int y1=0; y1<dimy; y1++)
      sumofsquares[x1]=(y1==0?0:sumofsquares[x1])+pow(theta[x1]-x[y1][x1],2);
  for (int x1=0; x1<dimx; x1++) {
    z[x1]=theta[x1]/stdErr(x1);
    }
  }
//---------------------------------------------------------------------------
double LogReg::stdErr(int idx) {
  return sqrt(sumofsquares[idx]/(dimy-1))/sqrt(dimy);
  }
//---------------------------------------------------------------------------
double LogReg::lowCI(int idx) {
  return theta[idx]-stdErr(idx)*CONF_LEVEL_95;
  }
//---------------------------------------------------------------------------
double LogReg::highCI(int idx) {
  return theta[idx]+stdErr(idx)*CONF_LEVEL_95;
  }
//---------------------------------------------------------------------------
bool LogReg::gradientDescent(int iterations, double minerror) {
  double cost1,cost2;

  clearArrays();
  cost1=calculateCost();
  for (int i1=0; i1<iterations; i1++) {
    calculateTheta();
    for (int x1=0; x1<dimx; x1++) {
      if (i1==0)
        theta[x1]=tmptheta[x1];
      else
        theta[x1]=theta[x1]-ALPHA*tmptheta[x1];
      }
    cost2=calculateCost();
    if (fabs(cost2-cost1)<minerror)
      return true;
    cost1=cost2;
    //cout << ">"<<Cost() << endl;
    //cout <<theta[0] << "\t" << theta[1]<< "\t" << theta[2] << "\t" << endl;
    //cout <<"**"<<tmptheta[0] << "\t" << tmptheta[1]<< "\t" << tmptheta[2] << "\t" << endl;

    }
  return false;
  }
//---------------------------------------------------------------------------
LogReg::~LogReg() {
  delete x;
  delete theta;
  delete tmptheta;
  delete sumofsquares;
  delete z;
  }
//---------------------------------------------------------------------------
