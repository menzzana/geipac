#include "logreg.h"
//---------------------------------------------------------------------------
LogisticRegression::~LogisticRegression() {  
  x.resize(0,0);
  variancecovariance.resize(0,0);
  y.resize(0);
  z.resize(0);
  oddsratio.resize(0);
  beta.resize(0);
  stderr.resize(0);
  }
//---------------------------------------------------------------------------
void LogisticRegression::clearArrays() {
  stderr.resize(x.cols());
  z.resize(x.cols());
  oddsratio.resize(x.cols());
  beta.resize(x.cols());
  beta.setConstant(x.cols(),0);
  variancecovariance.resize(0,0);
  }
//---------------------------------------------------------------------------
bool LogisticRegression::maximumLikelihoodRegression(int iterations, double minerror) {
  MatrixXd j,xt,jacobian;
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
      p(y1)=invLogit(p(y1));
      pnp=invOdds(p(y1));
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
    oddsratio(yx1)=exp(beta(yx1));
    }
  return iter1<iterations;
  }
//---------------------------------------------------------------------------
double LogisticRegression::invLogit(double p) {
  return exp(p)/(exp(p)+1);
  }
//---------------------------------------------------------------------------
double LogisticRegression::invOdds(double p) {
  return p*(1-p);
  }
//---------------------------------------------------------------------------
double LogisticRegression::lowCI(int idx) {
  return exp(beta(idx)-stderr(idx)*CONF_LEVEL_95);
  }
//---------------------------------------------------------------------------
double LogisticRegression::highCI(int idx) {
  return exp(beta(idx)+stderr(idx)*CONF_LEVEL_95);
  }
//---------------------------------------------------------------------------
double LogisticRegression::lowCI(double value, double error) {
  return value-error*CONF_LEVEL_95;
  }
//---------------------------------------------------------------------------
double LogisticRegression::highCI(double value, double error) {
  return value+error*CONF_LEVEL_95;
  }
//---------------------------------------------------------------------------
double LogisticRegression::calculateAPMValue(int idx1,int idx2,int idx3) {
  return exp(oddsratio(idx3)-oddsratio(idx1)-oddsratio(idx2)+1)/
         max(oddsratio(idx3),oddsratio(idx1)+oddsratio(idx2)-1);
  }
//------------------------------------------------------------------------------
double LogisticRegression::calculateRERI(int idx1,int idx2,int idx3) {
  return oddsratio(idx3)-oddsratio(idx1)-oddsratio(idx2)+1;
  }
//------------------------------------------------------------------------------
double LogisticRegression::APSEM(int idx1,int idx2,int idx3) {
  double ha1=-exp(beta(idx1)-beta(idx3));
  double ha2=-exp(beta(idx2)-beta(idx3));
  double ha3=(oddsratio(idx2)+oddsratio(idx1)-1)/oddsratio(idx3);
  return sqrt(pow(ha1,2)*variancecovariance(idx1,idx1)+
         pow(ha2,2)*variancecovariance(idx2,idx2)+
         pow(ha3,2)*variancecovariance(idx3,idx3)+
         2*ha1*ha2*variancecovariance(idx1,idx2)+
         2*ha1*ha3*variancecovariance(idx1,idx3)+
         2*ha2*ha3*variancecovariance(idx2,idx3));
  }
//------------------------------------------------------------------------------
