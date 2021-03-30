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

#include "datastore.h"
//------------------------------------------------------------------------------
using namespace GenEnvGen2I;
//------------------------------------------------------------------------------
DataStore::DataStore() {
  randomseed=CALC::rseed;
  appnegative=false;
  apcalculation=Proportion::DISEASE;
  threshold=THRESHOLD;
  iterations=ITERATIONS;
  permutations=0;
  naphenotype=0;
  cutoff=CUTOFF;
  model=Model::NONE;
  rawpermutation=false;
  markerid=nullptr;
  imarkerid=nullptr;
  individualid=nullptr;
  chromosome=nullptr;
  interactionfromfile=nullptr;
  covariate=nullptr;
  cutoff_mult=nullptr;
  cutoff_app=nullptr;
  permuted_mult=nullptr;
  permuted_app=nullptr;
  gender=nullptr;
  phenotype=nullptr;
  aphenotype=nullptr;
  genotype=nullptr;
  allele1=nullptr;
  allele2=nullptr;
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
  if (interactionfromfile!=nullptr)
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
  return (cutoff_app!=nullptr && cutoff_mult!=nullptr);
  }
//------------------------------------------------------------------------------
