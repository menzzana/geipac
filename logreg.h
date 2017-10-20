/*******************************************************************************
GEIPAC
Copyright (C) 2016  Henric Zazzi <henric@zazzi.se>

Geipac is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/

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
    double invLogit(double p);
    double invOdds(double p);
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
