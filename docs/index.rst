WHAT IS GEIPAC?
==============

This is the file README for the GEIPAC distribution

GEIPAC is a C++ parallel version of GEISA
GEISA is an update of JEIRA.
JEIRA is a Java implementation of GEIRA, a gene-environment interaction 
analysis tool written for R. It adds new capabilities, e.g. permutations
which was not present in the original implementation.

BUILDING (from source bundle)
=============================

**READ THIS CAREFULLY: you MUST install the required libraries!**

DEPENDENCIES
^^^^^^^^^^^^

In order to work, Geipac needs a couple of libraries.

==================== ===============================================================
C++                  Works well with the GNU compiler
cmake                version 2.8+
Boost                Version 1.36+ http://www.boost.org/
MPI                  If not found a serial version of parma will be compiled
==================== ===============================================================

Geipac has only been tested with the above mentioned versions, but may function with other versions as well.

In order to build just run::

  cd [binary dir]
  cmake [binary dir]/CMakeLists.txt
  make

USAGE
=====

geipac [OPTIONS]

OPTIONS
=======

:-n, --appnegative: Set this flag if negative APP values should
  be included in total permutation
:-a, --apcalculation <d/e/c>: Sets how the attributable proportion should be calculated.
  D calculates the proportion of the disease
  E calculates the proportion of the effect
  C calculates the corrected attributable proportion
  which takes negative values into consideration
  according to Hössjer et. al. Biostatistics,2014.
  Default: D
:-b,--basename <basename>: Specifies the base name of the binary
  input files (i.e. the name of the
  files without their file extensions:
  .bed, .bim, .fam).
:-c,--cutoff <n>: Specifies the minimum number of
  individuals in a group. The
  individuals are divided into groups
  (case/controls with and without the
  environmental factor). If any of these
  groups have a count below this value,
  no analysis will be performed on that
  marker. Default: 5
:-d,--model <type>: The model type to use (i.e. "dom" for
  dominant-dominant or rec for
  "recessive".
:-h,--help: Displays this help text.
:-i,--ifile <file>: Specifies the input interaction
  variable file. Default: null
:-f,--limitfile <file>: specifies a file containing
  significance limits for APp and MULT
  permutation calculations.
:-r,--iterations <iterations>: Sets the max number of iteration to
  perform when computing logistic
  regression (Default: 500)
:-t,--threshold <threshold>: Sets the min stable threshold when
  computing logistic regression
  (Default: 10E-3)
:-m,--markerfile <file>: Specifies a file containing
  interaction markers targeted for
  analysis.
:-o,--outputdir <path>: Specifies the directory where the
  output files will be stored. Default:
  None (Creates a result directory
  automatically)
:-p,--permutations <count>: Specifies the number of case/control
  permutations to perform. Default: 0
:-e,--permutationoutput <R=raw permutation output, T=total permutation output>:
  Sets if permutation rawdata should be
  printed to various files (Default: No)
  With Total permutation output, the first
  row shows the original data, whereas
  the other rows are the permutated data
:-s,--seed <value>: Specifies the seed used by the PRNG.
  Default: 123456789.

RECODE
======

If the presence of the risk allele is determined to be protective, a recode
is performed in such a way that what's considered a risk factor is reassessed. 
The recode is denoted in the output column "recode" with the following 
possible values:

0. No recode is performed.
1. The risk allele is considered to have a protective effect, and the risk 
   factor will be thought of as the absence of the risk allele. E.g. if 
   alleles A and T are present and A was initially considered the risk 
   allele, the absence of A will now be denoted the risk factor.
2. The interaction variable is inverted.
3. This is a combination of recode 1 and 2.

The column denoting the risk allele in the output will remain the same even 
after recoding. So in the case of a recode 1 or 3, the risk allele is in fact
considered protective.

UNIQUE FILE FORMAT
==================

Beside the standard binary input files, Geisa also contains other files
to interact with the data.

INTERACTION VARIABLE FILE
^^^^^^^^^^^^^^^^^^^^^^^^^

The interaction variable files contains individual IDs, environment variable and
covariates.
The first line of the file should depict the specific column data, and all
columns should be separated by TAB.
Individuals columns should be name INDID.
Environment variable should be name ENV
All other columns will be treated as covariate columns

Example.

INDID ENV COV1  TEST  HELLO
04D01801	0  1 0 1
 
First column is Individual ID, and 2nd is Environment.
COV1, TEST and HELLO are all covariates.
If no interaction variable file is present, the interaction will be calculated
from the genotype data.

LIMIT FILE
^^^^^^^^^^

The limit file contain only 2 columns.
The first line of the file should depict the specific column data, and all
columns should be separated by TAB.
The cutoff column for AP_pvalue should be named CUTOFF_APP
whereas the Multiplicative_interaction_term_pvalue cutoff column should
be name CUTOFF_MULT.
As many cutoff values as wanted can be added.

INTERACTION MARKER FILE
^^^^^^^^^^^^^^^^^^^^^^^

Should only contain one column with marker names.


COPYRIGHT
=========

GEIPAC is written by Henric Zazzi.
henric@zazzi.se


AVAILABILITY
============

The main web site for GEIPAC is https://bitbucket.org/menzzana/geipac
