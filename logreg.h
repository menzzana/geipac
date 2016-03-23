#ifndef LOGREG_H
#define LOGREG_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
class LogReg {
  public:
    static const double ALPHA=1;
    static const double CONF_LEVEL_95=1.96;
    double **x,**xo,*y,*theta,*tmptheta,*z,*sumofsquares;
    int dimx,dimy;

    LogReg();
    void createArrays(double *y,double **xo,int dimy,int dimx);
    void clearArrays();
    double matrixMultiply(double *x1,double *x2);
    double sigmoid(double z);
    void normalize();
    double calculateCost();
    void calculateTheta();
    void calculateZ();
    double stdErr(int idx);
    double lowCI(int idx);
    double highCI(int idx);
    bool gradientDescent(int iterations, double minerror);
    ~LogReg();
  };
//---------------------------------------------------------------------------
#endif // LOGREG_H
