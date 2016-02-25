#include "genenvgen2i.h"
//------------------------------------------------------------------------------
using namespace GenEnvGen2I;
//---------------------------------------------------------------------------
Analysis::Analysis() {
  param.randomseed=0;
  param.appnegative=false;
  param.apcalculation=NO_PROPORTION;
  param.threshold=THRESHOLD;
  param.iterations=ITERATIONS;
  param.permutations=0;
  param.cutoff=CUTOFF;
  param.model=NO_MODEL;
  param.permutation_output=PERMUTATION_TOTALDATA;
  }
//------------------------------------------------------------------------------
