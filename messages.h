#ifndef MESSAGES_H
#define MESSAGES_H
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
  const char *const SEED_OPTION[]={"seed,s","seed","Specifies the random seed used by the analysis [Default: 123456789]\n"};
  }
//------------------------------------------------------------------------------
// Error Messages
//------------------------------------------------------------------------------
namespace ERROR_TEXT {
  const char MPI_NOT_FOUND[]="Cannot initiate MPI";
  const char NO_MODEL_TYPE[]="No model type was defined. Must be either DOM or REC";
  const char NO_BED_FILE[]="Provided binary file is not a BED file";
  const char FILE_NOT_FOUND[]="File was not found: ";
  }
//------------------------------------------------------------------------------
namespace FILE_TEXT {
  const char RESULT[]="results.txt";
  const char MARKER_PERMUTATION_RESULT[]="marker_permutation_results.txt";
  const char TOTAL_PERMUTATION_RESULT[]="total_permutation_results.txt";
  const char TOTAL_PERMUTATIONS[]="total_permutations.txt";
  const char RESULT_PERMUTATION[]="results_permutation_%s.txt";
  }
//------------------------------------------------------------------------------
// Status Messages
//------------------------------------------------------------------------------
namespace STATUS_TEXT {
  const char IMARKER[]="Analyzing interaction with marker : %1$s";
  const char PERMUTATION_START[]="%tc Starting permutation iteration %d of %d";
  const char ORIGINAL_START[]="%tc Starting analysis of original results";
  const char WRITE_PERMUTATION[]="%tc: Writing permutation results";
  const char OUTPUT_READY[]="%tc: Ready.";
  const char COMPLETED_ITERATION[]="%tc: Completed iteration %d of %d";
  const char POSITIVE_CONVERGENCE[]="Yes";
  const char NEGATIVE_CONVERGENCE[]="No";
  const char FINISH[]="%tc Finished";
  }
//------------------------------------------------------------------------------
// Header Messages
//------------------------------------------------------------------------------
namespace HEADER_TEXT {
  const char RUN[]=        "Running JEIRA with the following parameters:";
  const char DATA_STORE[]= "Data store:               %s";
  const char FILE_BASE[]=  "File base:                %s";
  const char IFILE[]=      "Interaction file:         %s";
  const char IMARKERFILE[]="Interaction Marker file:  %s";
  const char LIMIT[]=      "Limit file:               %s";
  const char OUTPUT[]=     "Output directory:         %s";
  const char PERMUTATION[]="Permutations:             %d";
  const char SEED[]=       "Seed:                     %d";
  const char MODEL[]=      "Model Type:               %s";
  const char CUTOFF[]=     "Cutoff:                   %d";
  const char ITERATIONS[]= "LR Iterations:            %d";
  const char THRESHOLD[]=  "LR Threshold:             %f";
  const char APPNEG[]=     "Include APP negative:     %s";
  const char APCALC[]=     "AP Calculation method:    %s";
  }
//------------------------------------------------------------------------------
#endif // MESSAGES_H
