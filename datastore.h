#ifndef DATASTORE_H
#define DATASTORE_H

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "global.h"
//==============================================================================
namespace GenEnvGen2I {
  enum MODEL_TYPE { NO_MODEL, DOMINANT, RECESSIVE };
  enum PROPORTION_TYPE { NO_PROPORTION, PROPORTION_DISEASE, PROPORTION_EFFECT, PROPORTION_CORRECTED };
  
  static const double THRESHOLD=1E-3;
  static const int ITERATIONS=500;
  static const int CUTOFF=5;
  static const double MAX_P_VALUE=1;
  static const int ORIGINAL=0;
//------------------------------------------------------------------------------
  class DataStore {
    public:
      long randomseed;
      char apcalculation,model;
      int cutoff,iterations,permutations,naphenotype;
      int nimarkerid,nmarkerid,nlimit,nindividualid,ncovariate;
      double threshold;
      bool appnegative,rawpermutation;
      ostream *wres,*wperm,*wtotperm;
      
      string *markerid,*individualid,*chromosome;
      int *imarkerid,*gender,*interactionfromfile;
      int **phenotype,**genotype,**covariate,***aphenotype;
      double *cutoff_mult,*cutoff_app,*permuted_mult,*permuted_app;
      char *allele1,*allele2;
    
      DataStore();
      ~DataStore();
      void initialize();
      void permutePhenotypes();
      bool totalPermutations();
    };
//------------------------------------------------------------------------------
  }
//==============================================================================
#endif // DATASTORE_H
