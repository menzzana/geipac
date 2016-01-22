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
#include "boost/program_options.hpp"

//------------------------------------------------------------------------------
namespace prgm_opt=boost::program_options;
//------------------------------------------------------------------------------
void CleanUp(bool exitvalue) {
  exit(exitvalue);
  }
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
  prgm_opt::variables_map option_map;
  prgm_opt::options_description options("Options");

  try {
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
    prgm_opt::store(prgm_opt::parse_command_line(argc,argv,options),option_map);
    if (option_map.count(CMDOPTIONS::HELP_OPTION[1])) {
      printVersion();
      cout << options;
      CleanUp(EXIT_SUCCESS);
      }
    CleanUp(EXIT_SUCCESS);
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    CleanUp(EXIT_FAILURE);
    }
  }
//------------------------------------------------------------------------------

