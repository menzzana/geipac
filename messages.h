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

#ifndef MESSAGES_H
#define MESSAGES_H
//------------------------------------------------------------------------------
// strings used for options
//------------------------------------------------------------------------------
namespace CMDOPTIONS {
  const char *const HELP_OPTION[]={"help,h","help","Displays available options\n"};
  const char *const AP_OPTION[]={"apcalculation,a","apcalculation","Sets how the attributable proportion should be calculated. [D/E/C] Default: D\n"};
  const char *const BASE_OPTION[]={"basename,b","basename","Specify the base name of the plink input files.\n"};
  const char *const MODEL_OPTION[]={"model,d","model","The model type to use. [dom/rec]\n"};
  const char *const INTERACTION_OPTION[]={"interactionfile,i","interactionfile","Specifies the input interaction variable file.\n"};
  const char *const COVARIATE_OPTION[]={"covariatefile,v","covariatefile","Specifies the input covariate file.\n"};
  const char *const OUTPUT_OPTION[]={"outputdir,o","outputdir","Specifies the directory where the output files will be stored. If none is entered a result directory will be created automatically.\n"};
  const char *const PERMUTATION_OPTION[]={"permutations,p","permutations","Specifies the number of case/control permutations to perform. Default: 0\n"};
  const char *const LIMIT_OPTION[]={"limitfile,f","limitfile","specifies a file containing significance limits for APp and MULT permutation calculations.\n"};
  const char *const APPN_OPTION[]={"appnegative,n","appnegative","Set this flag if negative APP values should be included in total permutation.\n"};
  const char *const CUTOFF_OPTION[]={"cutoff,c","cutoff","Specifies the minimum number of individuals in a group. Default: 5\n"};
  const char *const ITERATION_OPTION[]={"iterations,r","iterations","Sets the max number of iteration to perform when computing logistic regression. Default: 500\n"};
  const char *const THRESHOLD_OPTION[]={"threshold,t","threshold","Sets the min stable threshold when computing logistic regression. Default: 10E-3\n"};
  const char *const MARKER_OPTION[]={"markerfile,m","markerfile","Specifies a file containing interaction markers targeted for analysis.\n"};
  const char *const RAWPERMUTATION_OPTION[]={"rawpermutation,w","rawpermutation","Sets if permutation rawdata should be written to result file.\n"};
  const char *const SEED_OPTION[]={"seed,s","seed","Specifies the random seed used by the analysis. Default: 123456789\n"};
  const char *const ALT_PHENOTYPE_OPTION[]={"phenotypes,y","phenotypes","Specifies a file for alternate mutiple phenotypes.\n"};
  const char *const LOG_OPTION[]={"log,l","Log output","Logs all output information in files\n"};
  }
//------------------------------------------------------------------------------
// Error Messages
//------------------------------------------------------------------------------
namespace ERROR_TEXT {
  const char NO_MODEL_TYPE[]="No model type was defined. Must be either DOM or REC";
  const char NO_BED_FILE[]="Provided binary file is not a BED file";
  const char NO_PLINK_FILES[]="Plink files data was not provided";
  const char UNKNOWN_INDIVIDUAL[]="Interaction variable/covariate individual is not present in dataset or missing individuals in interaction variable file or covariate file";
  const char MISSING_INDIVIDUAL[]="Individual %s not found.";
  const char NO_LIMITS[]="No limit columns could be found in file";
  const char NO_ALT_PHENOTYPE[]="No Alternate phenotype file could be found";
  const char FILE_NOT_FOUND[]="File was not found: %s";
  const char NO_FILE_LOADED[]="Wrong format of file: %s";
  const char MISSING_INTERACTION_MARKERS[]="Missing interaction marker(s)";
  const char MISSING_MARKER[]="Marker %s not found.";
  const char NO_OPTIONS[]="Try 'geipac --help' for more information.";
  }
//------------------------------------------------------------------------------
namespace FILE_TEXT {
  const char RESULT[]="/results.txt";
  const char MARKER_PERMUTATION_RESULT[]="/marker_permutation_results.txt";
  const char TOTAL_PERMUTATION_RESULT[]="/total_permutation_results.txt";
  const char OUTPUT_DIRECTORY[]="geipac_results_";
  const char LOG[]="/log.txt";
  }
//------------------------------------------------------------------------------
// Status Messages
//------------------------------------------------------------------------------
namespace STATUS_TEXT {
  const char IMARKER[]="Analyzing interaction with marker : %s";
  const char PERMUTATION_START[]=".";
  const char NO_PERMUTATION[]="No";
  const char ORIGINAL_START[]="#";
  const char WRITE_PERMUTATION[]="Writing permutation results";
  const char OUTPUT_READY[]="Ready.";
  const char COMPLETED_ITERATION[]="Completed iteration %d of %d";
  const char CONVERGENCE[]="Yes";
  const char NO_CONVERGENCE[]="No";
  const char FINISH[]="Finished";
  }
//------------------------------------------------------------------------------
// Header Messages
//------------------------------------------------------------------------------
namespace HEADER_TEXT {
  const char RUN[]=            "Running GEIPAC with the following parameters:";
  const char FILE_BASE[]=      "File base:                %s";
  const char INTERACTIONFILE[]="Interaction file:         %s";
  const char IMARKERFILE[]=    "Interaction Marker file:  %s";
  const char INTERACTION[]=    "Interaction:              %s";
  const char COVARIATE[]=      "Covariates:               %s";
  const char LIMIT[]=          "Limit file:               %s";
  const char ALTPHENOTYPE[]=   "Alternate Phenotype file: %s";
  const char OUTPUT[]=         "Output directory:         %s";
  const char PERMUTATION[]=    "Permutations:             %d";
  const char SEED[]=           "Seed:                     %d";
  const char MODEL[]=          "Model Type:               %s";
  const char CUTOFF[]=         "Cutoff:                   %d";
  const char ITERATIONS[]=     "LR Iterations:            %d";
  const char THRESHOLD[]=      "LR Threshold:             %f";
  const char APPNEG[]=         "Include APP negative:     %s";
  const char APCALC[]=         "AP Calculation method:    %s";
  const char PROCESSES[]=      "Processes:                %s";
  const char FROMGENEDATA[]="From genetic data";
  const char FROMVARFILE[]="From interaction variable file";
  const char NONE[]="None";
  const char YES[]="Yes";
  const char NO[]="No";
  }
//------------------------------------------------------------------------------
#endif // MESSAGES_H
