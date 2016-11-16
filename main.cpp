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
ofstream fpresult,fppermutation,fptotalpermutation;
//------------------------------------------------------------------------------
void CleanUp(bool exitvalue) {
  fpresult.close();
  fppermutation.close();
  fptotalpermutation.close();
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

  try {
    // Initialize
    imarker=NULL;
    ivariable=NULL;
    limit=NULL;
    plink=NULL;
    outputdir="";
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
      (CMDOPTIONS::RAWPERMUTATION_OPTION[0],CMDOPTIONS::RAWPERMUTATION_OPTION[2])
      (CMDOPTIONS::PERMUTATION_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::PERMUTATION_OPTION[2])
      (CMDOPTIONS::TOTALPERMUTATION_OPTION[0],CMDOPTIONS::TOTALPERMUTATION_OPTION[2])
      (CMDOPTIONS::SEED_OPTION[0],prgm_opt::value<double>()->required(),CMDOPTIONS::SEED_OPTION[2]);
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
    if (option_map.count(CMDOPTIONS::RAWPERMUTATION_OPTION[1]))
      myanalysis.param.rawpermutation=true;
    if (option_map.count(CMDOPTIONS::TOTALPERMUTATION_OPTION[1]))
      myanalysis.param.totalpermutation=true;
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
    // Print some information message.
    WRITELN(HEADER_TEXT::RUN);
    WRITELN_VALUE(HEADER_TEXT::FILE_BASE,global::getFileName(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()));
    WRITELN_VALUE(HEADER_TEXT::INTERACTIONFILE,(ivariable==NULL?"None":global::getFileName(option_map[CMDOPTIONS::INTERACTION_OPTION[1]].as<string>())));
    WRITELN_VALUE(HEADER_TEXT::IMARKERFILE,(imarker==NULL?"None":global::getFileName(option_map[CMDOPTIONS::MARKER_OPTION[1]].as<string>())));
    WRITELN_VALUE(HEADER_TEXT::INTERACTION,(imarker==NULL || ivariable==NULL?HEADER_TEXT::FROMGENEDATA:HEADER_TEXT::FROMVARFILE));
    WRITELN_VALUE(HEADER_TEXT::LIMIT,(limit==NULL?"None":global::getFileName(option_map[CMDOPTIONS::LIMIT_OPTION[1]].as<string>())));
    WRITELN_VALUE(HEADER_TEXT::OUTPUT,option_map[CMDOPTIONS::OUTPUT_OPTION[1]].as<string>());
    WRITELN_VALUE(HEADER_TEXT::PERMUTATION,myanalysis.param.permutations);
    WRITELN_VALUE(HEADER_TEXT::SEED,abs(myanalysis.param.randomseed));
    WRITELN_VALUE(HEADER_TEXT::MODEL,(myanalysis.param.model==GenEnvGen2I::DOMINANT?GenEnvGen2I::DOM_TEXT:GenEnvGen2I::REC_TEXT));
    WRITELN_VALUE(HEADER_TEXT::CUTOFF,myanalysis.param.cutoff);
    WRITELN_VALUE(HEADER_TEXT::ITERATIONS,myanalysis.param.iterations);
    WRITELN_VALUE(HEADER_TEXT::THRESHOLD,myanalysis.param.threshold);
    WRITELN_VALUE(HEADER_TEXT::APPNEG,(myanalysis.param.appnegative?"Yes":"No"));
    WRITELN_VALUE(HEADER_TEXT::APCALC,(myanalysis.param.apcalculation==GenEnvGen2I::DISEASE?GenEnvGen2I::DISEASE_TEXT:
      myanalysis.param.apcalculation==GenEnvGen2I::EFFECT?GenEnvGen2I::EFFECT_TEXT:GenEnvGen2I::CORRECTED_TEXT));
    // Transfer data to analysis class
    myanalysis.nindividualid=plink->fam->Length<FAMData>();
    myanalysis.nlimit=limit->Length<LimitData>();
    myanalysis.nmarkerid=plink->bim->Length<BIMData>();
    myanalysis.ncovariate=ivariable->ncovariate;
    myanalysis.initialize();
    if (imarker==NULL) {
      myanalysis.nimarkerid=myanalysis.nmarkerid;
      myanalysis.imarkerid=plink->bim->get<int>(&BIMData::index,myanalysis.nmarkerid,NULL);
      }
    else {
      myanalysis.nimarkerid=imarker->Length<IMarkerData>();
      myanalysis.imarkerid=imarker->get<int>(&IMarkerData::index,myanalysis.nimarkerid,NULL);
      }
    myanalysis.cutoff_app=limit->get<double>(&LimitData::cutoff_app,myanalysis.nlimit,NULL);
    myanalysis.cutoff_mult=limit->get<double>(&LimitData::cutoff_mult,myanalysis.nlimit,NULL);
    myanalysis.gender=plink->fam->get<int>(&FAMData::gender,myanalysis.nindividualid,NULL);
    plink->fam->get<int>(&FAMData::phenotype,myanalysis.nindividualid,myanalysis.phenotype[GenEnvGen2I::ORIGINAL]);
    myanalysis.individualid=plink->fam->get<string>(&FAMData::individualid,myanalysis.nindividualid,NULL);
    myanalysis.allele1=plink->bim->get<char>(&BIMData::allele1,myanalysis.nmarkerid,NULL);
    myanalysis.allele2=plink->bim->get<char>(&BIMData::allele2,myanalysis.nmarkerid,NULL);
    myanalysis.markerid=plink->bim->get<string>(&BIMData::markerid,myanalysis.nmarkerid,NULL);
    myanalysis.chromosome=plink->bim->get<string>(&BIMData::chromosome,myanalysis.nmarkerid,NULL);
    myanalysis.genotype=plink->getGenotypes(myanalysis.nindividualid,myanalysis.nmarkerid);
    myanalysis.covariate=ivariable->get(&IVariableData::covariate,myanalysis.nindividualid,myanalysis.ncovariate);
    if (ivariable->areInteractionsPresent() && imarker==NULL)
      myanalysis.interactionfromfile=ivariable->get<int>(&IVariableData::interaction,myanalysis.nindividualid,NULL);
    delete imarker;
    delete limit;
    delete plink->fam;
    delete plink->bim;
    delete plink;
    delete ivariable;
    // Output file headers
    fpresult.open((outputdir+FILE_TEXT::RESULT).c_str());
    myanalysis.param.wres=&fpresult;
    GenEnvGen2I::Analysis::printResults(*myanalysis.param.wres,RESULT_COLUMNS::TEXT,RESULT_COLUMNS::LENGTH_TEXT);
    GenEnvGen2I::Analysis::printResults(*myanalysis.param.wres,RESULT_COLUMNS::VALUES,RESULT_COLUMNS::LENGTH_VALUES);
    *myanalysis.param.wres<<endl;
    if (myanalysis.param.permutations>0) {
      fppermutation.open((outputdir+FILE_TEXT::MARKER_PERMUTATION_RESULT).c_str());
      myanalysis.param.wperm=&fppermutation;
      GenEnvGen2I::Analysis::printResults(*myanalysis.param.wperm,RESULT_COLUMNS::TEXT,RESULT_COLUMNS::PERM);
      GenEnvGen2I::Analysis::printResults(*myanalysis.param.wperm,RESULT_COLUMNS::VALUES,RESULT_COLUMNS::LENGTH_VALUES,RESULT_COLUMNS::PERMUTED_VALUE); 
      *myanalysis.param.wperm<<endl;
      if (myanalysis.param.totalpermutation) {
        fptotalpermutation.open((outputdir+FILE_TEXT::TOTAL_PERMUTATION_RESULT).c_str());
        myanalysis.param.wtotperm=&fptotalpermutation;
        GenEnvGen2I::Analysis::printResults(*myanalysis.param.wtotperm,RESULT_COLUMNS::TOTAL_PERMUTATIONS,RESULT_COLUMNS::LENGTH_TOTAL);
        *myanalysis.param.wtotperm<<endl;
        }
      }
    // Analysis
    myanalysis.createCovariateMatrix();
    if (myanalysis.param.permutations>0) {
      WRITELN(STATUS_TEXT::PERMUTE);
      myanalysis.permutePhenotypes();
      }
    for (int imarkeridx=0; imarkeridx<myanalysis.nimarkerid; imarkeridx++) {
      WRITELN_VALUE(STATUS_TEXT::IMARKER,myanalysis.markerid[myanalysis.imarkerid[imarkeridx]]);
      myanalysis.run(myanalysis.imarkerid[imarkeridx]);
      clog<<endl;
      }
    CleanUp(EXIT_SUCCESS);
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    CleanUp(EXIT_FAILURE);
    }
  }
//------------------------------------------------------------------------------
