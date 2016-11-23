#include "datastore.h"
//------------------------------------------------------------------------------
using namespace GenEnvGen2I;
//------------------------------------------------------------------------------
DataStore::DataStore() {
  randomseed=CALC::rseed;
  appnegative=false;
  apcalculation=NO_PROPORTION;
  threshold=THRESHOLD;
  iterations=ITERATIONS;
  permutations=0;
  cutoff=CUTOFF;
  model=NO_MODEL;
  rawpermutation=false;
  totalpermutation=false;
  markerid=NULL;
  imarkerid=NULL;
  individualid=NULL;
  chromosome=NULL;
  interactionfromfile=NULL;
  covariate=NULL;
  cutoff_mult=NULL;
  cutoff_app=NULL;
  permuted_mult=NULL;
  permuted_app=NULL;
  gender=NULL;
  phenotype=NULL;
  genotype=NULL;
  allele1=NULL;
  allele2=NULL;
  nimarkerid=0;
  nmarkerid=0;
  nindividualid=0;
  nlimit=0;
  ncovariate=0;
  }
//------------------------------------------------------------------------------
DataStore::~DataStore() {
  delete[] markerid;
  delete[] imarkerid;
  delete[] individualid;
  delete[] chromosome;
  if (interactionfromfile!=NULL)
    delete interactionfromfile;
  delete[] covariate;
  delete cutoff_mult;
  delete cutoff_app;
  delete permuted_mult;
  delete permuted_app;
  delete gender;
  delete[] phenotype;
  delete[] genotype;
  delete allele1;
  delete allele2;
  }
//------------------------------------------------------------------------------
void DataStore::initialize() {
  phenotype=global::make2DArray<int>(permutations+1,nindividualid);
  if (permutations>0) {
    permuted_mult=new double[permutations+1];
    permuted_app=new double[permutations+1];
    fill_n(permuted_mult,permutations+1,MAX_P_VALUE);
    fill_n(permuted_app,permutations+1,MAX_P_VALUE);
    }
  }
//------------------------------------------------------------------------------
void DataStore::permutePhenotypes() {
  CALC::sran1(randomseed);
  for (int permidx=1; permidx<=permutations; permidx++) {
    memcpy(phenotype[permidx],phenotype[permidx-1],nindividualid*sizeof(int));
    CALC::randomShuffle(phenotype[permidx],nindividualid);
    }
  }
//------------------------------------------------------------------------------