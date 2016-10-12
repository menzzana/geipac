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
    static constexpr double CONF_LEVEL_95=1.96;
    MatrixXd x;
    VectorXd y,beta,stderr,z;
    int firstx;

    void clearArrays();
    bool maximumLikelihoodRegression(int iterations, double minerror);
    double getMULTPropability(int idx);
    double lowCI(int idx);
    double highCI(int idx);
    ~LogisticRegression();
  };
//---------------------------------------------------------------------------
#endif // LOGREG_H
