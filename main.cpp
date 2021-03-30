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

#include "version_config.h"
#include "global.h"
#include "loader.h"
#include "datastore.h"
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
GenEnvGen2I::DataStore datastore;
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
  AltPhenotypeData *aphenotype;
  BEDData *plink;
  string outputdir;

  try {
    // Initialize
    imarker=nullptr;
    ivariable=nullptr;
    limit=nullptr;
    plink=nullptr;
    aphenotype=nullptr;
    outputdir="";
    // Program options
    prgm_opt::arg="[Value]";
    options.add_options()
      (CMDOPTIONS::HELP_OPTION[0],CMDOPTIONS::HELP_OPTION[2])
      (CMDOPTIONS::AP_OPTION[0],prgm_opt::value<char>()->required(),CMDOPTIONS::AP_OPTION[2])
      (CMDOPTIONS::BASE_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::BASE_OPTION[2])
      (CMDOPTIONS::MODEL_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::MODEL_OPTION[2])
      (CMDOPTIONS::INTERACTION_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::INTERACTION_OPTION[2])
      (CMDOPTIONS::OUTPUT_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::OUTPUT_OPTION[2])
      (CMDOPTIONS::PERMUTATION_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::PERMUTATION_OPTION[2])
      (CMDOPTIONS::LIMIT_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::LIMIT_OPTION[2])
      (CMDOPTIONS::APPN_OPTION[0],CMDOPTIONS::APPN_OPTION[2])
      (CMDOPTIONS::CUTOFF_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::CUTOFF_OPTION[2])
      (CMDOPTIONS::ITERATION_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::ITERATION_OPTION[2])
      (CMDOPTIONS::THRESHOLD_OPTION[0],prgm_opt::value<double>()->required(),CMDOPTIONS::THRESHOLD_OPTION[2])
      (CMDOPTIONS::MARKER_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::MARKER_OPTION[2])
      (CMDOPTIONS::RAWPERMUTATION_OPTION[0],CMDOPTIONS::RAWPERMUTATION_OPTION[2])
      (CMDOPTIONS::SEED_OPTION[0],prgm_opt::value<double>()->required(),CMDOPTIONS::SEED_OPTION[2])
      (CMDOPTIONS::ALT_PHENOTYPE_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::ALT_PHENOTYPE_OPTION[2]);
    prgm_opt::store(prgm_opt::parse_command_line(argc,argv,options),option_map);
    printVersion();
    if (option_map.count(CMDOPTIONS::HELP_OPTION[1])) {
      cout << options;
      CleanUp(EXIT_SUCCESS);
      }
    if (option_map.count(CMDOPTIONS::SEED_OPTION[1]))
      datastore.randomseed=-option_map[CMDOPTIONS::SEED_OPTION[1]].as<double>();
    if (option_map.count(CMDOPTIONS::CUTOFF_OPTION[1]))
      datastore.cutoff=option_map[CMDOPTIONS::CUTOFF_OPTION[1]].as<int>();
    if (option_map.count(CMDOPTIONS::ITERATION_OPTION[1]))
      datastore.iterations=option_map[CMDOPTIONS::ITERATION_OPTION[1]].as<int>();
    if (option_map.count(CMDOPTIONS::THRESHOLD_OPTION[1]))
      datastore.threshold=option_map[CMDOPTIONS::THRESHOLD_OPTION[1]].as<double>();
    if (option_map.count(CMDOPTIONS::PERMUTATION_OPTION[1]))
      datastore.permutations=option_map[CMDOPTIONS::PERMUTATION_OPTION[1]].as<int>();
    if (option_map.count(CMDOPTIONS::APPN_OPTION[1]))
      datastore.appnegative=true;
    if (option_map.count(CMDOPTIONS::MODEL_OPTION[1])) {
      string s1=boost::algorithm::to_lower_copy(option_map[CMDOPTIONS::MODEL_OPTION[1]].as<string>());
      datastore.model=s1.compare(GenEnvGen2I::DOM)==0?GenEnvGen2I::Model::DOMINANT:
        datastore.model=s1.compare(GenEnvGen2I::REC)==0?GenEnvGen2I::Model::RECESSIVE:GenEnvGen2I::Model::NONE;
      }
    if (option_map.count(CMDOPTIONS::RAWPERMUTATION_OPTION[1]))
      datastore.rawpermutation=true;
    if (option_map.count(CMDOPTIONS::AP_OPTION[1])) {
      char c1=tolower(option_map[CMDOPTIONS::AP_OPTION[1]].as<char>());
      datastore.apcalculation=(c1==GenEnvGen2I::EFFECT?GenEnvGen2I::Proportion::EFFECT:
        c1==GenEnvGen2I::CORRECTED?GenEnvGen2I::Proportion::CORRECTED:GenEnvGen2I::Proportion::DISEASE);
      }
    if (option_map.count(CMDOPTIONS::OUTPUT_OPTION[1]))
      outputdir=option_map[CMDOPTIONS::OUTPUT_OPTION[1]].as<string>();
    if (option_map.count(CMDOPTIONS::MARKER_OPTION[1]))
      imarker=IMarkerData::loadFile<IMarkerData>(option_map[CMDOPTIONS::MARKER_OPTION[1]].as<string>());
    if (option_map.count(CMDOPTIONS::INTERACTION_OPTION[1]))
      ivariable=IVariableData::loadFile<IVariableData>(option_map[CMDOPTIONS::INTERACTION_OPTION[1]].as<string>());
    if (option_map.count(CMDOPTIONS::LIMIT_OPTION[1])) {
      limit=LimitData::loadFile<LimitData>(option_map[CMDOPTIONS::LIMIT_OPTION[1]].as<string>());
      if (limit==nullptr)
        THROW_ERROR(ERROR_TEXT::NO_LIMITS);
      }
    if (option_map.count(CMDOPTIONS::BASE_OPTION[1])) {
      BIMData *bim;
      FAMData *fam;
      fam=nullptr;
      bim=nullptr;
      fam=FAMData::loadFile<FAMData>(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()+FAM_FILE);
      bim=BIMData::loadFile<BIMData>(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()+BIM_FILE);
      plink=BEDData::loadBinaryFile(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()+BED_FILE,fam,bim);
      }
    if (option_map.count(CMDOPTIONS::ALT_PHENOTYPE_OPTION[1])) {
      aphenotype=AltPhenotypeData::loadFile<AltPhenotypeData>(option_map[CMDOPTIONS::ALT_PHENOTYPE_OPTION[1]].as<string>());
      if (aphenotype==nullptr)
        THROW_ERROR(ERROR_TEXT::NO_ALT_PHENOTYPE);
      }      
    // Check received data
    if (plink==nullptr)
      THROW_ERROR(ERROR_TEXT::NO_PLINK_FILES);
    if (datastore.model==GenEnvGen2I::Model::NONE)
      THROW_ERROR(ERROR_TEXT::NO_MODEL_TYPE);
    if (!plink->bim->setInteractionMarkerIndex(imarker))
      THROW_ERROR(ERROR_TEXT::MISSING_INTERACTION_MARKERS);
    if (!ivariable->areAllIndividualPresent<IVariableData,FAMData>(plink->fam) ||
        !aphenotype->areAllIndividualPresent<AltPhenotypeData,FAMData>(plink->fam))
      THROW_ERROR(ERROR_TEXT::UNKNOWN_INDIVIDUAL);
    // set output
    outputdir=Loader::setOutputDirectory(outputdir);
    Loader::deleteResultFile(outputdir+FILE_TEXT::RESULT);
    Loader::deleteResultFile(outputdir+FILE_TEXT::MARKER_PERMUTATION_RESULT);
    Loader::deleteResultFile(outputdir+FILE_TEXT::TOTAL_PERMUTATION_RESULT);
    // Transfer data to analysis class
    datastore.nindividualid=plink->fam->Length<FAMData>();
    datastore.nlimit=limit->Length<LimitData>();
    datastore.nmarkerid=plink->bim->Length<BIMData>();
    datastore.ncovariate=ivariable==nullptr?0:ivariable->ncovariate;
    if (aphenotype!=nullptr)
      datastore.naphenotype=aphenotype->naphenotype;
    datastore.initialize();
    if (imarker==nullptr) {
      datastore.nimarkerid=datastore.nmarkerid;
      datastore.imarkerid=plink->bim->get<int>(&BIMData::index,datastore.nmarkerid,nullptr);
      }
    else {
      datastore.nimarkerid=imarker->Length<IMarkerData>();
      datastore.imarkerid=imarker->get<int>(&IMarkerData::index,datastore.nimarkerid,nullptr);
      }
    datastore.cutoff_app=limit->get<double>(&LimitData::cutoff_app,datastore.nlimit,nullptr);
    datastore.cutoff_mult=limit->get<double>(&LimitData::cutoff_mult,datastore.nlimit,nullptr);
    datastore.gender=plink->fam->get<int>(&FAMData::gender,datastore.nindividualid,nullptr);
    plink->fam->get<int>(&FAMData::phenotype,datastore.nindividualid,datastore.phenotype[GenEnvGen2I::ORIGINAL]);
    datastore.individualid=plink->fam->get<string>(&FAMData::individualid,datastore.nindividualid,nullptr);
    datastore.allele1=plink->bim->get<char>(&BIMData::allele1,datastore.nmarkerid,nullptr);
    datastore.allele2=plink->bim->get<char>(&BIMData::allele2,datastore.nmarkerid,nullptr);
    datastore.markerid=plink->bim->get<string>(&BIMData::markerid,datastore.nmarkerid,nullptr);
    datastore.chromosome=plink->bim->get<string>(&BIMData::chromosome,datastore.nmarkerid,nullptr);
    datastore.genotype=plink->getGenotypes(datastore.nindividualid,datastore.nmarkerid);
    datastore.covariate=ivariable->get<int>(&IVariableData::covariate,datastore.nindividualid,datastore.ncovariate,nullptr);
    if (aphenotype!=nullptr)
      aphenotype->get<int>(&AltPhenotypeData::aphenotype,datastore.nindividualid,datastore.naphenotype,datastore.aphenotype[GenEnvGen2I::ORIGINAL]);
    if (ivariable->areInteractionsPresent() && imarker==nullptr)
      datastore.interactionfromfile=ivariable->get<int>(&IVariableData::interaction,datastore.nindividualid,nullptr);
    // Print some information message.
    WRITELN(HEADER_TEXT::RUN);
    WRITELN_VALUE(HEADER_TEXT::FILE_BASE,global::getFileName(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()));
    WRITELN_VALUE(HEADER_TEXT::INTERACTIONFILE,(ivariable==nullptr?"None":global::getFileName(option_map[CMDOPTIONS::INTERACTION_OPTION[1]].as<string>())));
    WRITELN_VALUE(HEADER_TEXT::IMARKERFILE,(imarker==nullptr?"None":global::getFileName(option_map[CMDOPTIONS::MARKER_OPTION[1]].as<string>())));
    WRITELN_VALUE(HEADER_TEXT::INTERACTION,(datastore.interactionfromfile==nullptr?HEADER_TEXT::FROMGENEDATA:HEADER_TEXT::FROMVARFILE));
    WRITELN_VALUE(HEADER_TEXT::LIMIT,(limit==nullptr?"None":global::getFileName(option_map[CMDOPTIONS::LIMIT_OPTION[1]].as<string>())));
    WRITELN_VALUE(HEADER_TEXT::OUTPUT,outputdir);
    WRITELN_VALUE(HEADER_TEXT::PERMUTATION,datastore.permutations);
    WRITELN_VALUE(HEADER_TEXT::SEED,abs(datastore.randomseed));
    WRITELN_VALUE(HEADER_TEXT::MODEL,(datastore.model==GenEnvGen2I::Model::DOMINANT?GenEnvGen2I::DOM_TEXT:GenEnvGen2I::REC_TEXT));
    WRITELN_VALUE(HEADER_TEXT::CUTOFF,datastore.cutoff);
    WRITELN_VALUE(HEADER_TEXT::ITERATIONS,datastore.iterations);
    WRITELN_VALUE(HEADER_TEXT::THRESHOLD,datastore.threshold);
    WRITELN_VALUE(HEADER_TEXT::APPNEG,(datastore.appnegative?"Yes":"No"));
    WRITELN_VALUE(HEADER_TEXT::APCALC,(datastore.apcalculation==GenEnvGen2I::Proportion::DISEASE?GenEnvGen2I::DISEASE_TEXT:
      datastore.apcalculation==GenEnvGen2I::Proportion::EFFECT?GenEnvGen2I::EFFECT_TEXT:GenEnvGen2I::CORRECTED_TEXT));
    WRITELN_VALUE(HEADER_TEXT::ALTPHENOTYPE,(aphenotype==nullptr?"None":global::getFileName(option_map[CMDOPTIONS::ALT_PHENOTYPE_OPTION[1]].as<string>())));
    // Delete pointers structures for import of data
    imarker->Delete<IMarkerData>();
    limit->Delete<LimitData>();
    plink->fam->Delete<FAMData>();
    plink->bim->Delete<BIMData>();
    plink->Delete<BEDData>();
    ivariable->Delete<IVariableData>();
    aphenotype->Delete<AltPhenotypeData>();
    // Output file headers
    fpresult.open((outputdir+FILE_TEXT::RESULT).c_str());
    datastore.wres=&fpresult;
    GenEnvGen2I::Analysis::printResults(*datastore.wres,RESULT_COLUMNS::TEXT,RESULT_COLUMNS::LENGTH_TEXT);
    GenEnvGen2I::Analysis::printResults(*datastore.wres,RESULT_COLUMNS::VALUES,RESULT_COLUMNS::LENGTH_VALUES);
    *datastore.wres<<endl;
    if (datastore.permutations>0) {
      fppermutation.open((outputdir+FILE_TEXT::MARKER_PERMUTATION_RESULT).c_str());
      datastore.wperm=&fppermutation;
      GenEnvGen2I::Analysis::printResults(*datastore.wperm,RESULT_COLUMNS::TEXT,RESULT_COLUMNS::PERM);
      GenEnvGen2I::Analysis::printResults(*datastore.wperm,RESULT_COLUMNS::VALUES,RESULT_COLUMNS::LENGTH_VALUES,RESULT_COLUMNS::PERMUTED_VALUE);
      *datastore.wperm<<endl;
      if (datastore.totalPermutations()) {
        fptotalpermutation.open((outputdir+FILE_TEXT::TOTAL_PERMUTATION_RESULT).c_str());
        datastore.wtotperm=&fptotalpermutation;
        GenEnvGen2I::Analysis::printResults(*datastore.wtotperm,RESULT_COLUMNS::TOTAL_PERMUTATIONS,RESULT_COLUMNS::LENGTH_TOTAL);
        *datastore.wtotperm<<endl;
        }
      }
    // Analysis
    if (datastore.permutations>0)
      datastore.permutePhenotypes();
    #pragma omp parallel
      {
      #ifndef SERIAL
      if (omp_get_thread_num()==0)
        WRITELN_VALUE(HEADER_TEXT::PROCESSES,omp_get_num_threads());
      #endif
      GenEnvGen2I::Analysis *myanalysis;
      myanalysis=new GenEnvGen2I::Analysis(&datastore);
      for (int imarkeridx=0; imarkeridx<datastore.nimarkerid; imarkeridx++) {
      #ifndef SERIAL
      if (omp_get_thread_num()==0)
      #endif
        WRITELN_VALUE(STATUS_TEXT::IMARKER,datastore.markerid[datastore.imarkerid[imarkeridx]]);
      myanalysis->run(datastore.imarkerid[imarkeridx]);
      #ifndef SERIAL
      if (omp_get_thread_num()==0)
      #endif
        clog<<endl;
        }
      delete myanalysis;
      }
    GenEnvGen2I::Analysis::printTotalPermutation(datastore);    
    CleanUp(EXIT_SUCCESS);
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    CleanUp(EXIT_FAILURE);
    }
  }
//------------------------------------------------------------------------------
