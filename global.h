/*******************************************************************************
GEIPAC
Copyright (C) 2016  Henric Zazzi <hzazzi@kth.se>

This program is free software: you can redistribute it and/or modify
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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "messages.h"

#ifndef SERIAL
  #include "mpi.h"
#endif

#define THROW_ERROR(text) throw runtime_error(global::to_string(text))
#define THROW_ERROR_VALUE(text,value) throw runtime_error(global::to_string(boost::format(text) % value))
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
// global constants and functions used
//------------------------------------------------------------------------------
namespace global {
  static const int MPIROOT=0;
  template<typename T> string to_string(T value) {  //lexical_cast does funny things with double
    stringstream s1;

    s1 << value;
    return s1.str();
    }
  }
//------------------------------------------------------------------------------
namespace CALC {
  #define IA 16807
  #define IM 2147483647
  #define AM (1.0/IM)
  #define IQ 127773
  #define IR 2836
  #define NTAB 32
  #define NDIV (1+(IM-1)/NTAB)
  #define EPS 1.2e-7
  #define RNMX (1.0-EPS)
  static long rseed=-123456789;

  void sran1(long rseed);
  double ran1();
  }
//------------------------------------------------------------------------------
// strings used for options
//------------------------------------------------------------------------------
namespace CMDOPTIONS {
  const char *const HELP_OPTION[]={"help,h","help","Displays available commands\n"};
  const char *const APPN_OPTION[]={"appnegative,n","appnegative","if negative APP values should be included in total permutation calculations\n"};
  const char *const AP_OPTION[]={"apcalculation,a","apcalculation","Sets how the attributable proportion should be calculated. [D/E/C]\n"};
  const char *const BASE_OPTION[]={"basename,b","basename","Specify the base name of the binary input files\n"};
  const char *const CUTOFF_OPTION[]={"cutoff,c","cutoff","Specifies the minimum number of individuals in a group [Default: 5]\n"};
  const char *const MODEL_OPTION[]={"model,d","model","The model type to use [dom/rec]\n"};
  const char *const INTERACTION_OPTION[]={"interactionfile,i","interactionfile","Specifies the input interaction variable file\n"};
  const char *const LIMIT_OPTION[]={"limitfile,f","limitfile","specifies a file containing significance limits for APp and MULT permutation calculations\n"};
  const char *const ITERATION_OPTION[]={"iterations,r","iterations","Sets the max number of iteration to perform when computing logistic regression [Default: 500]\n"};
  const char *const THRESHOLD_OPTION[]={"threshold,t","threshold","Sets the min stable threshold when computing logistic regression [Default: 10E-3]\n"};
  const char *const MARKER_OPTION[]={"markerfile,m","markerfile","Specifies a file containing interaction markers targeted for analysis.\n"};
  const char *const OUTPUT_OPTION[]={"outputdir,o","outputdir","Specifies the directory where the output files will be stored. Default: None (Creates a result directory automatically)\n"};
  const char *const PERMUTATION_OPTION[]={"permutations,p","permutations","Specifies the number of case/control permutations to perform. Default: 0\n"};
  const char *const PERMUTATIONOUTPUT_OPTION[]={"permutationoutput,e","permutationouput","Sets if permutation rawdata should be printed to various files [R/T]\n"};
  const char *const SEED_OPTION[]={"seed,s","seed","Specifies the random seed used by the analysis (Default: 123456789]\n"};
  }
//------------------------------------------------------------------------------
#endif // GLOBAL_H
