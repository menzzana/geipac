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

#ifndef GENENVGEN2I_H
#define GENENVGEN2I_H

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "global.h"
#include "datastore.h"
#include "logreg.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
//==============================================================================
namespace GenEnvGen2I {
  enum class Phenotype {
    UNKNOWN, UNAFFECTED,AFFECTED
    };
  enum class Gender {
    UNKNOWN,MALE,FEMALE
    };
  enum class Allele {
    CONTROL_PRIMARY,CONTROL_SECONDARY,CASE_PRIMARY,CASE_SECONDARY,OFFSET=1,NA=0
    };
  enum class Ap {
    A0B0,A1B0,A0B1,A1B1,LENGTH
    };
  enum class Apm {
    A1,B1,A1B1,LENGTH
    };
  enum class Lr {
    BETA0,A1B0,A0B1,A1B1,LENGTH
    };
  enum class Lrm {
    BETA0,A1,B1,A1B1,LENGTH
    };
  enum class Zygosity {
    HOMOZYGOTE_PRIMARY, UNKNOWN, HETEROZYGOTE, HOMOZYGOTE_SECONDARY
    };
  enum class RiskFactor {
    NO,YES,NA=-1
    };
  enum class Interaction {
    NO,YES,NA=-1
    };
  #define CHROMOSOME_X "X"
  #define DELIMITER "\t"
  static const char DISEASE='d';
  static const char EFFECT='e';
  static const char CORRECTED='c';
  static const char DISEASE_TEXT[]="Disease";
  static const char EFFECT_TEXT[]="Effect";
  static const char CORRECTED_TEXT[]="Corrected";
  static const int N_RISK_MATRIX=8;
  static const char REC[]="rec";
  static const char DOM[]="dom";
  static const char REC_TEXT[]="Recessive";
  static const char DOM_TEXT[]="Dominant";
  static const char TOTAL_PERMUTATION[]="Total";
//------------------------------------------------------------------------------
  class Analysis {
    public:
      DataStore *data;
      double **covariate1,**covariate2;
      int *interaction,*imarkinteraction;
      RiskFactor *riskfactors;
      
      Analysis(DataStore *datastore);
      ~Analysis();
      void setInteraction(int interactivemarkeridx);
      void run(int interactivemarkeridx);
      static void printTotalPermutation(DataStore &data1);
      void analyzeData(int markeridx, int *phenotypex,string *results_text, double *results_value);
      void alleleSummaryCount(int *alleles,int markeridx,int * phenotypex);
      bool validIndividualData(int individualidx,int markeridx,int *phenotypex);
      bool validGeneticData(int individualidx,int markeridx,int *phenotypex);
      char calculateRiskAllele(int markeridx, int *phenotypex, string *results);
      void calculateRiskFactors(int markeridx,int *phenotypex,char riskallele,int recode);
      bool isDominantOrXMale(int individualidx,int markeridx);
      void calculateRiskMatrix(int *phenotypex,int *riskmatrix);
      bool belowCutOff(int *riskmatrix);
      void setCleanData(int markeridx,int *y, double **x, VectorXd &desty, MatrixXd &destx, int dimx, bool a0b0);
      void swapInteractions();
//------------------------------------------------------------------------------
      template<typename T> static void printResults(ostream &stream, T *results, int length) {	
        for (int i1=0; i1<length; i1++)
          stream << results[i1] << DELIMITER;
        }
//------------------------------------------------------------------------------
      template<typename T> static void printResults(ostream &stream, T *results, int length, bool const *print) {	
        for (int i1=0; i1<length; i1++)
          if (print[i1])
            stream << results[i1] << DELIMITER;
        }
//-----------------------------------------------------------------------------
    };
  }
//==============================================================================
namespace RESULT_COLUMNS {
  const char *const TEXT[]={
    "Interaction_marker","Chr_test_marker","Test_marker", "Permutation",
    "Test_marker_minor_allele","Test_marker_major_allele", "Test_marker_risk_allele"
    };

  enum INDEX_TEXT {
    INTERACTION, CHR, SNP, PERM, MINOR, MAJOR, RISK,
    LENGTH_TEXT
    };
//------------------------------------------------------------------------------
  const char *const VALUES[]={
    "ORa_double_exposure", "ORa_double_exposure_lower_limit","ORa_double_exposure_higher_limit", "ORa_test_marker",
    "ORa_test_marker_lower_limit", "ORa_test_marker_higher_limit", "ORa_risk_factor", "ORa_risk_factor_lower_limit",
    "ORa_risk_factor_higher_limit", "AP", "AP_L", "AP_H", "AP_pvalue", "Stable_additive_logistic_regression",
    "Multiplicative_interaction_term_pvalue", "ORm_interaction", "ORm_interaction_L", "ORm_interaction_H",
    "ORm_testmarker", "ORm_testmarker_L", "ORm_testmarker_H", "ORm_riskfactor", "ORm_riskfactor_L", "ORm_riskfactor_H",
    "APM", "APM_L", "APM_H", "APM_pvalue", "Stable_multiplicative_logistic_regression", "No_controls_test_0_risk_0",
    "No_cases_test_0_risk_0", "No_controls_test_0_risk_1", "No_cases_test_0_risk_1", "No_controls_test_1_risk_0",
    "No_cases_test_1_risk_0", "No_controls_test_1_risk_1", "No_cases_test_1_risk_1", "Recode_code"
    };

  enum INDEX_VALUES {
    ORII, ORIIL, ORIIH, ORIO, ORIOL, ORIOH, OROI, OROIL, OROIH, AP, APL, APH,	APP,
    STABLELRA, MULT, ORMII, ORMIIL, ORMIIH, ORMIO, ORMIOL, ORMIOH, ORMOI, ORMOIL,
    ORMOIH, APM, APML, APMH, APMP, STABLELRM, IND00_0, IND00_1, IND01_0, IND01_1,
    IND10_0, IND10_1, IND11_0, IND11_1, RECODE,
    LENGTH_VALUES
    };
  bool const PERMUTED_VALUE[] {
    true, true, true, true, true, true, true, true, true, true, false, false, true,
    true, true, false, false, false, false, false, false, false, false,
    false,true, false, false, false, true, false, false, false, false,
    false, false, false, false, false
    };
//------------------------------------------------------------------------------
  enum INDEX_TOTAL_PERMUTATIONS {
    TPERM,TSIGNAPP,TAPP,TSIGNMULT,TMULT,
    LENGTH_TOTAL
    };
    
  const char *const TOTAL_PERMUTATIONS[]={
    "Permutation","Significance Limit", "APP_permutation_pvalue",
    "Significance Limit", "MULT_permutation_pvalue"
    };
//------------------------------------------------------------------------------
  }
//==============================================================================
#endif // GENENVGEN2I_H
