#ifndef GENENVGEN2I_H
#define GENENVGEN2I_H

#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
//------------------------------------------------------------------------------
namespace GenEnvGen2I {

  enum PROPORTION { NO_PROPORTION, PROPORTION_DISEASE, PROPORTION_EFFECT, PROPORTION_CORRECTED };
  enum MODEL { NO_MODEL, DOMINANT, RECESSIVE };
  enum PERMUTATION_OUTPUT { PERMUTATION_RAWDATA, PERMUTATION_TOTALDATA };
  static const char DISEASE='d';
  static const char EFFECT='e';
  static const char CORRECTED='c';
  static const double THRESHOLD=1E-3;
  static const int ITERATIONS=500;
  static const int CUTOFF=5;
  static const char REC[]="rec";
  static const char DOM[]="dom";
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

      Analysis();
    };
  }
//------------------------------------------------------------------------------
#endif // GENENVGEN2I_H
