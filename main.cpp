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

#include "version_config.h"
#include "global.h"
#include "loader.h"
#include "logreg.h"
#include "genenvgen2i.h"
#include "boost/program_options.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctype.h>
#include <string>

//------------------------------------------------------------------------------
namespace prgm_opt=boost::program_options;
//------------------------------------------------------------------------------
GenEnvGen2I::Analysis myanalysis;
//------------------------------------------------------------------------------
void CleanUp(bool exitvalue) {
  exit(exitvalue);
  }
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
  prgm_opt::variables_map option_map;
  prgm_opt::options_description options("Options");
  IMarkerData *imarker;
  IVariableData *ivariable;
  LimitData *limit;
  BEDData *plink;
  string outputdir;
  int mpirank,mpisize;

  try {
    // Initialize
    imarker=NULL;
    ivariable=NULL;
    limit=NULL;
    plink=NULL;
    outputdir="";
    #ifndef SERIAL
      if (MPI_Init(&argc,&argv)!=MPI_SUCCESS)
        THROW_ERROR(ERROR_TEXT::MPI_NOT_FOUND);
      MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
      MPI_Comm_size(MPI_COMM_WORLD,&mpisize);
    #else
      mpirank=global::MPIROOT;
      mpisize=1;
    #endif
    // Program options
    prgm_opt::arg="[Value]";
    options.add_options()
      (CMDOPTIONS::HELP_OPTION[0],CMDOPTIONS::HELP_OPTION[2])
      (CMDOPTIONS::APPN_OPTION[0],CMDOPTIONS::APPN_OPTION[2])
      (CMDOPTIONS::AP_OPTION[0],prgm_opt::value<char>()->required(),CMDOPTIONS::AP_OPTION[2])
      (CMDOPTIONS::BASE_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::BASE_OPTION[2])
      (CMDOPTIONS::CUTOFF_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::CUTOFF_OPTION[2])
      (CMDOPTIONS::MODEL_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::MODEL_OPTION[2])
      (CMDOPTIONS::INTERACTION_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::INTERACTION_OPTION[2])
      (CMDOPTIONS::LIMIT_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::LIMIT_OPTION[2])
      (CMDOPTIONS::ITERATION_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::ITERATION_OPTION[2])
      (CMDOPTIONS::THRESHOLD_OPTION[0],prgm_opt::value<double>()->required(),CMDOPTIONS::THRESHOLD_OPTION[2])
      (CMDOPTIONS::MARKER_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::MARKER_OPTION[2])
      (CMDOPTIONS::OUTPUT_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::OUTPUT_OPTION[2])
      (CMDOPTIONS::PERMUTATION_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::PERMUTATION_OPTION[2])
      (CMDOPTIONS::PERMUTATIONOUTPUT_OPTION[0],prgm_opt::value<char>()->required(),CMDOPTIONS::PERMUTATIONOUTPUT_OPTION[2])
      (CMDOPTIONS::SEED_OPTION[0],prgm_opt::value<double>()->required(),CMDOPTIONS::SEED_OPTION[2]);
    if (mpirank==global::MPIROOT) {
      prgm_opt::store(prgm_opt::parse_command_line(argc,argv,options),option_map);
      printVersion();
      if (option_map.count(CMDOPTIONS::HELP_OPTION[1])) {
        cout << options;
        CleanUp(EXIT_SUCCESS);
        }
      if (option_map.count(CMDOPTIONS::SEED_OPTION[1]))
        myanalysis.param.randomseed=-option_map[CMDOPTIONS::SEED_OPTION[1]].as<double>();
      if (option_map.count(CMDOPTIONS::CUTOFF_OPTION[1]))
        myanalysis.param.cutoff=option_map[CMDOPTIONS::CUTOFF_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::ITERATION_OPTION[1]))
        myanalysis.param.iterations=option_map[CMDOPTIONS::ITERATION_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::THRESHOLD_OPTION[1]))
        myanalysis.param.threshold=option_map[CMDOPTIONS::THRESHOLD_OPTION[1]].as<double>();
      if (option_map.count(CMDOPTIONS::PERMUTATION_OPTION[1]))
        myanalysis.param.permutations=option_map[CMDOPTIONS::PERMUTATION_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::APPN_OPTION[1]))
        myanalysis.param.appnegative=true;
      if (option_map.count(CMDOPTIONS::MODEL_OPTION[1])) {
        string s1=boost::algorithm::to_lower_copy(option_map[CMDOPTIONS::MODEL_OPTION[1]].as<string>());
        myanalysis.param.model=s1.compare(GenEnvGen2I::DOM)==0?GenEnvGen2I::DOMINANT:
          myanalysis.param.model=s1.compare(GenEnvGen2I::REC)==0?GenEnvGen2I::RECESSIVE:0;
        }
      if (option_map.count(CMDOPTIONS::PERMUTATIONOUTPUT_OPTION[1])) {
        char c1=tolower(option_map[CMDOPTIONS::PERMUTATIONOUTPUT_OPTION[1]].as<char>());
        myanalysis.param.permutation_output=(c1==GenEnvGen2I::RAWDATA?GenEnvGen2I::PERMUTATION_RAWDATA:
          c1==GenEnvGen2I::TOTALDATA?GenEnvGen2I::PERMUTATION_TOTALDATA:0);
        }
      if (option_map.count(CMDOPTIONS::AP_OPTION[1])) {
        char c1=tolower(option_map[CMDOPTIONS::AP_OPTION[1]].as<char>());
        myanalysis.param.apcalculation=(c1==GenEnvGen2I::DISEASE?GenEnvGen2I::PROPORTION_DISEASE:
          c1==GenEnvGen2I::EFFECT?GenEnvGen2I::PROPORTION_EFFECT:
          c1==GenEnvGen2I::CORRECTED?GenEnvGen2I::PROPORTION_CORRECTED:0);
        }
      if (option_map.count(CMDOPTIONS::OUTPUT_OPTION[1]))
        outputdir=option_map[CMDOPTIONS::OUTPUT_OPTION[1]].as<string>();
      if (option_map.count(CMDOPTIONS::MARKER_OPTION[1]))
        imarker=IMarkerData::loadFile<IMarkerData>(option_map[CMDOPTIONS::MARKER_OPTION[1]].as<string>());
      if (option_map.count(CMDOPTIONS::INTERACTION_OPTION[1]))
        ivariable=IVariableData::loadFile<IVariableData>(option_map[CMDOPTIONS::INTERACTION_OPTION[1]].as<string>());
      if (option_map.count(CMDOPTIONS::LIMIT_OPTION[1])) {
        limit=LimitData::loadFile<LimitData>(option_map[CMDOPTIONS::LIMIT_OPTION[1]].as<string>());
        if (limit==NULL)
          THROW_ERROR(ERROR_TEXT::NO_LIMITS);
        }
      if (option_map.count(CMDOPTIONS::BASE_OPTION[1])) {
        BIMData *bim;
        FAMData *fam;
        fam=NULL;
        bim=NULL;
        fam=FAMData::loadFile<FAMData>(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()+FAM_FILE);
        bim=BIMData::loadFile<BIMData>(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()+BIM_FILE);
        plink=BEDData::loadBinaryFile(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()+BED_FILE,fam,bim);
        }
      // Check received data
      if (plink==NULL)
        THROW_ERROR(ERROR_TEXT::NO_PLINK_FILES);
      if (myanalysis.param.model==GenEnvGen2I::NO_MODEL)
        THROW_ERROR(ERROR_TEXT::NO_MODEL_TYPE);
      if (!plink->bim->setInteractionMarkerIndex(imarker))
        THROW_ERROR(ERROR_TEXT::MISSING_INTERACTION_MARKERS);
      if (!ivariable->areAllIndividualPresent(plink->fam))
        THROW_ERROR(ERROR_TEXT::UNKNOWN_INDIVIDUAL);
      // set output
      outputdir=Loader::setOutputDirectory(outputdir);
      Loader::deleteResultFile(FILE_TEXT::RESULT);
      Loader::deleteResultFile(FILE_TEXT::MARKER_PERMUTATION_RESULT);
      Loader::deleteResultFile(FILE_TEXT::TOTAL_PERMUTATION_RESULT);
      Loader::deleteResultFile(FILE_TEXT::TOTAL_PERMUTATIONS);
      for (int i1=1; Loader::deleteResultFile(global::to_string(boost::format(FILE_TEXT::RESULT_PERMUTATION) % i1)); i1++);
      // Print some information message.
      WRITE(HEADER_TEXT::RUN);
      WRITE_VALUE(HEADER_TEXT::FILE_BASE,global::getFileName(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()));
      WRITE_VALUE(HEADER_TEXT::INTERACTIONFILE,(ivariable==NULL?"None":global::getFileName(option_map[CMDOPTIONS::INTERACTION_OPTION[1]].as<string>())));
      WRITE_VALUE(HEADER_TEXT::IMARKERFILE,(imarker==NULL?"None":global::getFileName(option_map[CMDOPTIONS::MARKER_OPTION[1]].as<string>())));
      WRITE_VALUE(HEADER_TEXT::INTERACTION,(imarker==NULL || ivariable==NULL?HEADER_TEXT::FROMGENEDATA:HEADER_TEXT::FROMVARFILE));
      WRITE_VALUE(HEADER_TEXT::LIMIT,(limit==NULL?"None":global::getFileName(option_map[CMDOPTIONS::LIMIT_OPTION[1]].as<string>())));
      WRITE_VALUE(HEADER_TEXT::OUTPUT,option_map[CMDOPTIONS::OUTPUT_OPTION[1]].as<string>());
      WRITE_VALUE(HEADER_TEXT::PERMUTATION,myanalysis.param.permutations);
      WRITE_VALUE(HEADER_TEXT::SEED,abs(myanalysis.param.randomseed));
      WRITE_VALUE(HEADER_TEXT::MODEL,(myanalysis.param.model==GenEnvGen2I::DOMINANT?GenEnvGen2I::DOM_TEXT:GenEnvGen2I::REC_TEXT));
      WRITE_VALUE(HEADER_TEXT::CUTOFF,myanalysis.param.cutoff);
      WRITE_VALUE(HEADER_TEXT::ITERATIONS,myanalysis.param.iterations);
      WRITE_VALUE(HEADER_TEXT::THRESHOLD,myanalysis.param.threshold);
      WRITE_VALUE(HEADER_TEXT::APPNEG,(myanalysis.param.appnegative?"Yes":"No"));
      WRITE_VALUE(HEADER_TEXT::APCALC,(myanalysis.param.apcalculation==GenEnvGen2I::DISEASE?GenEnvGen2I::DISEASE_TEXT:
        myanalysis.param.apcalculation==GenEnvGen2I::EFFECT?GenEnvGen2I::EFFECT_TEXT:GenEnvGen2I::CORRECTED_TEXT));
      // Transfer data to analysis class
      myanalysis.nindividualid=plink->fam->Length<FAMData>();
      myanalysis.nlimit=limit->Length<LimitData>();
      myanalysis.nmarkerid=plink->bim->Length<BIMData>();
      if (imarker==NULL) {
        myanalysis.nimarkerid=myanalysis.nmarkerid;
        myanalysis.imarkerid=plink->bim->get(&BIMData::index,myanalysis.nmarkerid);
        }
      else {
        myanalysis.nimarkerid=imarker->Length<IMarkerData>();
        myanalysis.imarkerid=imarker->get(&IMarkerData::index,myanalysis.nimarkerid);
        }
      myanalysis.cutoff_app=limit->get(&LimitData::cutoff_app,myanalysis.nlimit);
      myanalysis.cutoff_mult=limit->get(&LimitData::cutoff_mult,myanalysis.nlimit);
      myanalysis.gender=plink->fam->get(&FAMData::gender,myanalysis.nindividualid);
      myanalysis.phenotype=plink->fam->get(&FAMData::phenotype,myanalysis.nindividualid);
      myanalysis.individualid=plink->fam->get(&FAMData::individualid,myanalysis.nindividualid);
      myanalysis.allele1=plink->bim->get(&BIMData::allele1,myanalysis.nmarkerid);
      myanalysis.allele2=plink->bim->get(&BIMData::allele2,myanalysis.nmarkerid);
      myanalysis.markerid=plink->bim->get(&BIMData::markerid,myanalysis.nmarkerid);
      myanalysis.chromosome=plink->bim->get(&BIMData::chromosome,myanalysis.nmarkerid);
      myanalysis.genotype=plink->getGenotypes(myanalysis.nindividualid,myanalysis.nmarkerid);
      if (ivariable->areInteractionsPresent() && imarker==NULL)
        myanalysis.interactionfromfile=ivariable->get(&IVariableData::interaction,myanalysis.nindividualid);
      myanalysis.ncovariate=ivariable->ncovariate;
      myanalysis.covariate=ivariable->get(&IVariableData::covariate,myanalysis.nindividualid,myanalysis.ncovariate);
      delete imarker;
      delete limit;
      delete plink->fam;
      delete plink->bim;
      delete plink;
      delete ivariable;
      }
    // Analysis
    myanalysis.initialize();
    GenEnvGen2I::Analysis::printResults(cout,NULL);
    for (int imarkeridx=0; imarkeridx<myanalysis.nimarkerid; imarkeridx++) {
      WRITE_VALUE(STATUS_TEXT::IMARKER,myanalysis.markerid[myanalysis.imarkerid[imarkeridx]]);
      myanalysis.run(myanalysis.imarkerid[imarkeridx]);
      }
    CleanUp(EXIT_SUCCESS);
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    CleanUp(EXIT_FAILURE);
    }
  }
//------------------------------------------------------------------------------

/*
ifstream fp;
int i1,f1,f2;
double **x,*y;
string fstr;

y=new double[100];
x=new double *[100];
x[0]=new double[200];
for (int i1=1; i1<100; i1++)
  x[i1]=&x[0][i1*2];
fp.open("LGOctave/ex2data1.txt");
i1=0;
while (getline(fp,fstr)) {
  f1=fstr.find_first_of(",");
  f2=fstr.find_last_of(",");
  y[i1]=strtod(fstr.substr(f2+1,f2).c_str(),NULL);
  x[i1][0]=strtod(fstr.substr(0,f1).c_str(),NULL);
  x[i1][1]=strtod(fstr.substr(f1+1,f2-f1-1).c_str(),NULL);
  i1++;
  }
fp.close();
LogReg lr=LogReg();
lr.createArrays(y,x,100,2);
lr.clearArrays();
lr.Normalize();
cout << lr.Cost() << endl;
cout << "------------------------------------------------"<<endl;
lr.GradientDescent(400,1E-5);
cout << "------------------------------------------------"<<endl;
cout << lr.Cost() << endl;
cout << lr.theta[0] << "\t" << lr.theta[1]<< "\t" << lr.theta[2] << "\t" << endl;
exit(EXIT_SUCCESS);
*/
