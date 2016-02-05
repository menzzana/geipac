#ifndef LOGREG_H
#define LOGREG_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
//---------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
class LogReg {
  private:
    static const double ALPHA=1;
  public:
    double **x,**xo,*y,*theta,*tmptheta;
    int dimx,dimy;

    LogReg();
    void createArrays(double *y,double **xo,int dimy,int dimx);
    void clearArrays();
    double MatrixMultiply(double *x1,double *x2);
    double sigmoid(double z);
    void Normalize();
    double Cost();
    void CalculateTheta();
    bool GradientDescent(int iterations, double minerror);
    ~LogReg();
  };
//---------------------------------------------------------------------------
#endif // LOGREG_H
