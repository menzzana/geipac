#ifndef GENENVGEN2I_H
#define GENENVGEN2I_H

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <ctype.h>
#include "global.h"
#include "logreg.h"
//------------------------------------------------------------------------------
namespace GenEnvGen2I {

  enum PROPORTION_TYPE { NO_PROPORTION, PROPORTION_DISEASE, PROPORTION_EFFECT, PROPORTION_CORRECTED };
  enum MODEL_TYPE { NO_MODEL, DOMINANT, RECESSIVE };
  enum PERMUTATION_OUTPUT_TYPE { PERMUTATION_RAWDATA, PERMUTATION_TOTALDATA };
  enum PHENOTYPE_TYPE {PHENOTYPE_UNKNOWN, PHENOTYPE_UNAFFECTED,PHENOTYPE_AFFECTED};
  enum GENDER_TYPE {GENDER_UNKNOWN,GENDER_MALE,GENDER_FEMALE};
  enum INDEX_TYPE {INDEX_CONTROL_PRIMARY,INDEX_CONTROL_SECONDARY,INDEX_CASE_PRIMARY,INDEX_CASE_SECONDARY,INDEX_POSITION_OFFSET=1,INDEX_NA=0};
  enum MATRIX_TYPE {LRTHETA0,MATRIX_INDEX_A0B0,MATRIX_INDEX_A1B0,MATRIX_INDEX_A0B1,MATRIX_INDEX_A1B1,MATRIX_INDEX_COV2};
  enum MATRIX_MULT_TYPE {LRTHETA0_M,MATRIX_INDEX_A1m,MATRIX_INDEX_B1m,MATRIX_INDEX_A1mB1m,MATRIX_INDEX_COV1};
  enum ZYGOSITY_TYPE {HOMOZYGOTE_PRIMARY, ZYGOTE_UNKNOWN, HETEROZYGOTE, HOMOZYGOTE_SECONDARY};
  enum RISKFACTOR_TYPE {NO_RISK,RISK,NA_RISK=-1};
  enum INTERACTION_TYPE {NO_INTERACTION,INTERACTION,NA_INTERACTION=-1};
  #define CHROMOSOME_X "X"
  static const char DISEASE='d';
  static const char EFFECT='e';
  static const char CORRECTED='c';
  static const char DISEASE_TEXT[]="Disease";
  static const char EFFECT_TEXT[]="Effect";
  static const char CORRECTED_TEXT[]="Corrected";
  static const double THRESHOLD=1E-3;
  static const int ITERATIONS=500;
  static const int CUTOFF=5;
  static const char REC[]="rec";
  static const char DOM[]="dom";
  static const char REC_TEXT[]="Recessive";
  static const char DOM_TEXT[]="Dominant";
  static const char RAWDATA='r';
  static const char TOTALDATA='t';

//------------------------------------------------------------------------------
  class Analysis {
    public:
      struct Param {
        long randomseed;
        char apcalculation,model,permutation_output;
        int cutoff,iterations,permutations;
        double threshold;
        bool appnegative;
        } param;

      string *markerid,*individualid,*chromosome;
      int *interaction,*interactionfromfile,*imarkinteraction,**covariate;
      double *cutoff_mult,*cutoff_app;
      int *gender,*phenotype,*permphenotype,**genotype;
      char *allele1,*allele2;
      int *imarkerid,nimarkerid,nmarkerid,nlimit,nindividualid,ncovariate,*riskfactors;
      double **covdata1,**covdata2,**cleancovdata;
      double *rephenotype,*cleanrephenotype;
      LogisticRegression logreg1,logreg2;

      Analysis();
      ~Analysis();
      void setInteraction(int interactivemarkeridx);
      void initialize();
      void run(int interactivemarkeridx);
      void alleleSummaryCount(int *alleles,int markeridx);
      bool validIndividualData(int individualidx,int markeridx);
      bool validGeneticData(int individualidx,int markeridx);
      char calculateRiskAllele(int markeridx, string *results);
      void calculateRiskFactors(int markeridx,char riskallele,int recode);
      bool isDominantOrXMale(int individualidx,int markeridx);
      bool calculateRiskMatrix(string *results);
      int cleanData(int markeridx,double *y, double **x, int dimx);
      void shufflePhenotype();
      void swapInteractions();
      static void printResults(ostream &stream,string *results);
    };
  }
//------------------------------------------------------------------------------
namespace RESULT_COLUMNS {
  const char *const RESULT_COLUMN_TEXT[]={"perm","Interaction_marker","Chr_test_marker","Test_marker",
    "ORa_double_exposure", "ORa_double_exposure_lower_limit","ORa_double_exposure_higher_limit", "ORa_test_marker",
    "ORa_test_marker_lower_limit", "ORa_test_marker_higher_limit", "ORa_risk_factor", "ORa_risk_factor_lower_limit",
    "ORa_risk_factor_higher_limit", "AP", "AP_L", "AP_H", "AP_pvalue", "Stable_additive_logistic_regression",
    "Multiplicative_interaction_term_pvalue", "ORm_interaction", "ORm_interaction_L", "ORm_interaction_H",
    "ORm_testmarker", "ORm_testmarker_L", "ORm_testmarker_H", "ORm_riskfactor", "ORm_riskfactor_L", "ORm_riskfactor_H",
    "APM", "APM_L", "APM_H", "APM_pvalue", "Stable_multiplicative_logistic_regression", "No_controls_test_0_risk_0",
    "No_cases_test_0_risk_0", "No_controls_test_0_risk_1", "No_cases_test_0_risk_1", "No_controls_test_1_risk_0",
    "No_cases_test_1_risk_0",  "No_controls_test_1_risk_1", "No_cases_test_1_risk_1", "Test_marker_minor_allele",
    "Test_marker_major_allele", "Test_marker_risk_allele", "recode_code", "Temporary_threshold"
    };

  enum RESULT_COLUMN_IDX {PERM, INTERACTION, CHR, SNP, ORII, ORIIL, ORIIH, ORIO, ORIOL, ORIOH, OROI,
    OROIL, OROIH, AP, APL, APH,	APP, STABLELRA, MULT, ORMII, ORMIIL, ORMIIH, ORMIO, ORMIOL, ORMIOH,
    ORMOI, ORMOIL, ORMOIH, APM, APML, APMH, APMP, STABLELRM, IND00_0, IND00_1, IND01_0, IND01_1,
    IND10_0, IND10_1,	IND11_0, IND11_1,	MINOR, MAJOR, RISK, RECODE, THRESHOLD
    };
  static const int LENGTH_RESULTS=46;
  }
//------------------------------------------------------------------------------
#endif // GENENVGEN2I_H
