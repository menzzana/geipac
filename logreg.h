#ifndef LOGREG_H
#define LOGREG_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
class LogisticRegression {
  public:
    static constexpr double ALPHA=1;
    static constexpr double CONF_LEVEL_95=1.96;
    double **x,*y,*theta,*tmptheta,*z,*sumofsquares;
    int dimx,dimy;

    LogisticRegression();
    void createArrays(double *y,double **x,int dimx);
    void clearArrays();
    double matrixMultiply(double *x1,double *x2);
    double sigmoid(double z);
    double logit(double z);
    double getMULTPropability(int idx);
    double calculateCost();
    void calculateTheta();
    void calculateZ();
    double stdErr(int idx);
    double lowCI(int idx);
    double highCI(int idx);
    bool gradientDescent(int iterations, double minerror);
    ~LogisticRegression();
  };
//---------------------------------------------------------------------------
#endif // LOGREG_H
