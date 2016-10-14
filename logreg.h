#ifndef LOGREG_H
#define LOGREG_H

#include <stdio.h>
#include "global.h"
#include <iostream>
#include <string.h>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/LU>
//------------------------------------------------------------------------------
using namespace std;
using namespace Eigen;
//------------------------------------------------------------------------------
// Logistic regression class
//------------------------------------------------------------------------------
class LogisticRegression {
  public:
    MatrixXd x,variancecovariance;
    VectorXd y,beta,stderr,z,oddsratio;

    void clearArrays();
    bool maximumLikelihoodRegression(int iterations, double minerror);
    double lowCI(int idx);
    double highCI(int idx);
    double lowCI(double value, double error);
    double highCI(double value, double error);
    double calculateAPMValue(int idx1,int idx2,int idx3);
    double calculateRERI(int idx1,int idx2,int idx3);
    double APSEM(int idx1,int idx2,int idx3);
    ~LogisticRegression();

  private:
    static constexpr double CONF_LEVEL_95=1.96;
  };
//---------------------------------------------------------------------------
#endif // LOGREG_H
