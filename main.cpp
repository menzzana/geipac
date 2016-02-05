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
#include "logreg.h"
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
void CleanUp(bool exitvalue) {
  exit(exitvalue);
  }
//------------------------------------------------------------------------------
int main(int argc, char **argv) {



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
      (CMDOPTIONS::CUTOFF_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::CUTOFF_OPTION[2])
      (CMDOPTIONS::MODEL_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::MODEL_OPTION[2])
      (CMDOPTIONS::INTERACTION_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::INTERACTION_OPTION[2])
      (CMDOPTIONS::LIMIT_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::LIMIT_OPTION[2])
      (CMDOPTIONS::ITERATION_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::ITERATION_OPTION[2])
      (CMDOPTIONS::THRESHOLD_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::THRESHOLD_OPTION[2])
      (CMDOPTIONS::MARKER_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::MARKER_OPTION[2])
      (CMDOPTIONS::OUTPUT_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::OUTPUT_OPTION[2])
      (CMDOPTIONS::PERMUTATION_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::PERMUTATION_OPTION[2])
      (CMDOPTIONS::PERMUTATIONOUTPUT_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::PERMUTATIONOUTPUT_OPTION[2])
      (CMDOPTIONS::SEED_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::SEED_OPTION[2]);
    if (mpirank==global::MPIROOT) {
      prgm_opt::store(prgm_opt::parse_command_line(argc,argv,options),option_map);
      if (option_map.count(CMDOPTIONS::HELP_OPTION[1])) {
        printVersion();
        cout << options;
        CleanUp(EXIT_SUCCESS);
        }

      }

    CleanUp(EXIT_SUCCESS);
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    CleanUp(EXIT_FAILURE);
    }
  }
//------------------------------------------------------------------------------

