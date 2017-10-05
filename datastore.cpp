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
  naphenotype=0;
  cutoff=CUTOFF;
  model=NO_MODEL;
  rawpermutation=false;
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
  aphenotype=NULL;
  genotype=NULL;
  allele1=NULL;
  allele2=NULL;
  nimarkerid=0;
  nmarkerid=0;
  nindividualid=0;
  nlimit=0;
  ncovariate=0;
  naphenotype=0;
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
  delete[] aphenotype;
  delete allele1;
  delete allele2;
  }
//------------------------------------------------------------------------------
void DataStore::initialize() {
  phenotype=global::make2DArray<int>(permutations+1,nindividualid);
  if (naphenotype>0)
    aphenotype=global::make3DArray<int>(permutations+1,nindividualid,naphenotype);
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
    if (naphenotype==0)
      continue;
    memcpy(aphenotype[permidx][0],aphenotype[permidx-1][0],naphenotype*nindividualid*sizeof(int));
    for (int apidx=0; apidx<naphenotype; apidx++)
      CALC::randomShuffle(aphenotype[permidx][apidx],nindividualid);
    }
  }
//------------------------------------------------------------------------------
bool DataStore::totalPermutations() {
  return (cutoff_app!=NULL && cutoff_mult!=NULL);
  }
//------------------------------------------------------------------------------
