#include "logreg.h"
#include "global.h"
//---------------------------------------------------------------------------
LogisticRegression::LogisticRegression() {
  theta=tmptheta=z=sumofsquares=NULL;
  dimx=dimy=0;
  }
//---------------------------------------------------------------------------
LogisticRegression::~LogisticRegression() {
  delete theta;
  delete tmptheta;
  delete sumofsquares;
  delete z;
  }
//---------------------------------------------------------------------------
void LogisticRegression::createArrays(double *y,double **x,int dimx) {
  this->dimx=dimx;
  this->y=y;
  this->x=x;
  theta=new double[this->dimx];
  tmptheta=new double[this->dimx];
  z=new double[this->dimx];
  sumofsquares=new double[this->dimx];
  }
//---------------------------------------------------------------------------
void LogisticRegression::clearArrays() {
  fill_n(theta,dimx,0);
  fill_n(z,dimx,0);
  fill_n(sumofsquares,dimx,0);
  }
//---------------------------------------------------------------------------
double LogisticRegression::matrixMultiply(double *n1,double *n2) {
  double res;

  res=0;
  for (int i1=0; i1<dimx; i1++)
    res+=n1[i1]*n2[i1];
  return res;
  }
//---------------------------------------------------------------------------
double LogisticRegression::sigmoid(double z) {
  return 1/(1+exp(-z));
  }
//---------------------------------------------------------------------------
double LogisticRegression::logit(double z) {
  return exp(z)/(exp(z)+1);
  }
//---------------------------------------------------------------------------
double LogisticRegression::getMULTPropability(int idx) {
  return 1-CALC::Chi2(sumofsquares[idx]/pow(z[idx],2),dimy);
  }
//------------------------------------------------------------------------------
  double LogisticRegression::calculateCost() {
  double J,h0;

  J=0;
  for (int y1=0; y1<dimy; y1++) {
    h0=sigmoid(matrixMultiply(x[y1],theta));
    J=J-y[y1]*log(h0)-(1-y[y1])*log(1-h0);
    }
  return J/dimy;
  }
//---------------------------------------------------------------------------
void LogisticRegression::calculateTheta() {
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
void LogisticRegression::calculateZ() {
  for (int x1=0; x1<dimx; x1++)
    for (int y1=0; y1<dimy; y1++)
      sumofsquares[x1]+=pow(theta[x1]-x[y1][x1],2);
  for (int x1=0; x1<dimx; x1++) {
    z[x1]=theta[x1]/stdErr(x1);
    }
  }
//---------------------------------------------------------------------------
double LogisticRegression::stdErr(int idx) {
  return sqrt(sumofsquares[idx]/(dimy-1))/sqrt(dimy);
  }
//---------------------------------------------------------------------------
double LogisticRegression::lowCI(int idx) {
  return theta[idx]-stdErr(idx)*CONF_LEVEL_95;
  }
//---------------------------------------------------------------------------
double LogisticRegression::highCI(int idx) {
  return theta[idx]+stdErr(idx)*CONF_LEVEL_95;
  }
//---------------------------------------------------------------------------
bool LogisticRegression::gradientDescent(int iterations, double minerror) {

  /*
  for (int i1=0; i1<iterations; i1++) {
    for (int x1=0; x1<dimx; x1++) {
      tmpbeta[x1]=beta[x1];
      p[x1]=logit(beta[x1]);
      }

    }

*/
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
    }
  return false;
  }
//---------------------------------------------------------------------------

