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
GenEnvGen2I::Analysis *myanalysis;
//------------------------------------------------------------------------------
void CleanUp(bool exitvalue) {
  delete myanalysis;
  exit(exitvalue);
  }
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
  prgm_opt::variables_map option_map;
  prgm_opt::options_description options("Options");
  int mpirank,mpisize;

  try {
    #ifndef SERIAL
      if (MPI_Init(&argc,&argv)!=MPI_SUCCESS)
        THROW_ERROR(ERROR_TEXT::MPI_NOT_FOUND);
      MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
      MPI_Comm_size(MPI_COMM_WORLD,&mpisize);
    #else
      mpirank=global::MPIROOT;
      mpisize=1;
    #endif
    prgm_opt::arg="[Value]";
    options.add_options()
      (CMDOPTIONS::HELP_OPTION[0],CMDOPTIONS::HELP_OPTION[2])
      (CMDOPTIONS::APPN_OPTION[0],CMDOPTIONS::APPN_OPTION[2])
      (CMDOPTIONS::AP_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::AP_OPTION[2])
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
    myanalysis=new GenEnvGen2I::Analysis();
    if (mpirank==global::MPIROOT) {
      prgm_opt::store(prgm_opt::parse_command_line(argc,argv,options),option_map);
      printVersion();
      if (option_map.count(CMDOPTIONS::HELP_OPTION[1])) {
        cout << options;
        CleanUp(EXIT_SUCCESS);
        }
      if (option_map.count(CMDOPTIONS::SEED_OPTION[1]))
        myanalysis->param.randomseed=-option_map[CMDOPTIONS::SEED_OPTION[1]].as<double>();
      if (option_map.count(CMDOPTIONS::CUTOFF_OPTION[1]))
        myanalysis->param.cutoff=option_map[CMDOPTIONS::CUTOFF_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::ITERATION_OPTION[1]))
        myanalysis->param.iterations=option_map[CMDOPTIONS::ITERATION_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::THRESHOLD_OPTION[1]))
        myanalysis->param.threshold=option_map[CMDOPTIONS::THRESHOLD_OPTION[1]].as<double>();
      if (option_map.count(CMDOPTIONS::PERMUTATION_OPTION[1]))
        myanalysis->param.permutations=option_map[CMDOPTIONS::PERMUTATION_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::APPN_OPTION[1]))
        myanalysis->param.appnegative=true;
      if (option_map.count(CMDOPTIONS::MODEL_OPTION[1])) {
        string s1=boost::algorithm::to_lower_copy(option_map[CMDOPTIONS::MODEL_OPTION[1]].as<string>());
        myanalysis->param.model=s1.compare(GenEnvGen2I::DOM)==0?GenEnvGen2I::DOMINANT:
          myanalysis->param.model=s1.compare(GenEnvGen2I::REC)==0?GenEnvGen2I::RECESSIVE:0;
        }
      if (option_map.count(CMDOPTIONS::PERMUTATIONOUTPUT_OPTION[1])) {
        char c1=tolower(option_map[CMDOPTIONS::PERMUTATIONOUTPUT_OPTION[1]].as<char>());
        myanalysis->param.permutation_output=(c1==GenEnvGen2I::RAWDATA?GenEnvGen2I::PERMUTATION_RAWDATA:
          c1==GenEnvGen2I::TOTALDATA?GenEnvGen2I::PERMUTATION_TOTALDATA:0);
        }
      if (option_map.count(CMDOPTIONS::AP_OPTION[1])) {
        char c1=tolower(option_map[CMDOPTIONS::AP_OPTION[1]].as<char>());
        myanalysis->param.apcalculation=(c1==GenEnvGen2I::DISEASE?GenEnvGen2I::PROPORTION_DISEASE:
          c1==GenEnvGen2I::EFFECT?GenEnvGen2I::PROPORTION_EFFECT:
          c1==GenEnvGen2I::CORRECTED?GenEnvGen2I::PROPORTION_CORRECTED:0);
        }
      if (option_map.count(CMDOPTIONS::OUTPUT_OPTION[1])) {
        string outputdir=option_map[CMDOPTIONS::OUTPUT_OPTION[1]].as<string>();

        }



      if (option_map.count(CMDOPTIONS::MARKER_OPTION[1])) {
        IMarkerData *imark;
        imark=NULL;
        imark=IMarkerData::loadFile(imark,option_map[CMDOPTIONS::MARKER_OPTION[1]].as<string>());

        IMarkerData *ld;
        for (ld=imark; ld!=NULL; ld=ld->Next)
          cout << ld->markerid << endl;

        }
      if (option_map.count(CMDOPTIONS::INTERACTION_OPTION[1])) {
        IVariableData *ivariable;
        ivariable=NULL;
        ivariable=IVariableData::loadFile(ivariable,option_map[CMDOPTIONS::INTERACTION_OPTION[1]].as<string>());

        }
      if (option_map.count(CMDOPTIONS::LIMIT_OPTION[1])) {
        LimitData *limit;
        limit=NULL;
        limit=LimitData::loadFile(limit,option_map[CMDOPTIONS::LIMIT_OPTION[1]].as<string>());

        }
      if (option_map.count(CMDOPTIONS::BASE_OPTION[1])) {
        FAMData *fam;
        BIMData *bim;
        BEDData *plink;
        fam=NULL;
        bim=NULL;
        fam=FAMData::loadFile(fam,option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()+FAM_FILE);
        bim=BIMData::loadFile(bim,option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()+BIM_FILE);
        plink=BEDData::loadBinaryFile(option_map[CMDOPTIONS::BASE_OPTION[1]].as<string>()+BED_FILE,fam,bim);

        }

      }
    if (myanalysis->param.model==GenEnvGen2I::NO_MODEL)
      THROW_ERROR(ERROR_TEXT::NO_MODEL_TYPE);

		// Deleting previous result files
		global::deleteResultFile(FILE_TEXT::RESULT);
		global::deleteResultFile(FILE_TEXT::MARKER_PERMUTATION_RESULT);
		global::deleteResultFile(FILE_TEXT::TOTAL_PERMUTATION_RESULT);
		global::deleteResultFile(FILE_TEXT::TOTAL_PERMUTATIONS);
		for (int i1=1; global::deleteResultFile(global::to_string(boost::format(FILE_TEXT::RESULT_PERMUTATION) % i1)); i1++);

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
