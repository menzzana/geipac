#include "logreg.h"
//---------------------------------------------------------------------------
LogisticRegression::~LogisticRegression() {  
  x.resize(0,0);
  y.resize(0);
  z.resize(0);
  beta.resize(0);
  stderr.resize(0);
  }
//---------------------------------------------------------------------------
void LogisticRegression::clearArrays() {
  stderr.resize(x.cols());
  z.resize(x.cols());
  beta.resize(x.cols());
  beta.setConstant(x.cols(),0);
  }
//---------------------------------------------------------------------------
bool LogisticRegression::maximumLikelihoodRegression(int iterations, double minerror) {
  MatrixXd j,xt,jacobian,variancecovariance;
  VectorXd oldbeta,p,s;
  double pnp;
  int iter1;

  clearArrays();
  j.resize(x.cols(),x.rows());
  xt=x.transpose();
  for (iter1=0; iter1<iterations; iter1++) {
    oldbeta=beta;
    p=x*beta;
    for (int y1=0; y1<x.rows(); y1++) {
      p(y1)=exp(p(y1))/(exp(p(y1))+1);
      pnp=p(y1)*(1-p(y1));
      p(y1)=y(y1)-p(y1);
      for (int x1=0; x1<x.cols(); x1++)
        j(x1,y1)=x(y1,x1)*pnp;
      }
    s=xt*p;
    jacobian=j*x;
    FullPivLU<MatrixXd> lu(jacobian);
    variancecovariance=lu.inverse();
    beta=beta+variancecovariance*s;
    if ((beta-oldbeta).array().abs().sum()<minerror)
      break;
    }
  for (int yx1=0; yx1<stderr.rows(); yx1++) {
    stderr(yx1)=sqrt(variancecovariance(yx1,yx1));
    z(yx1)=beta(yx1)/stderr(yx1);
    }
  return iter1==iterations;
  }
//------------------------------------------------------------------------------
double LogisticRegression::getMULTPropability(int idx) {
  return 1-CALC::Chi2(pow(z(idx),2),1);
  }
//---------------------------------------------------------------------------
double LogisticRegression::lowCI(int idx) {
  return beta(idx)-stderr(idx)*CONF_LEVEL_95;
  }
//---------------------------------------------------------------------------
double LogisticRegression::highCI(int idx) {
  return beta(idx)+stderr(idx)*CONF_LEVEL_95;
  }
//---------------------------------------------------------------------------
