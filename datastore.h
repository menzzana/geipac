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
  enum class Model { NONE,DOMINANT,RECESSIVE };
  enum class Proportion { DISEASE,EFFECT,CORRECTED };
  
  static const double THRESHOLD=1E-3;
  static const int ITERATIONS=500;
  static const int CUTOFF=5;
  static const double MAX_P_VALUE=1;
  static const int ORIGINAL=0;
//------------------------------------------------------------------------------
  class DataStore {
    public:
      long randomseed;
      Proportion apcalculation;
      Model model;
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
