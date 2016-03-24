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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "messages.h"
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/local_time/local_time.hpp>

#ifndef SERIAL
  #include "mpi.h"
#endif

#define THROW_ERROR(text) throw runtime_error(global::to_string(text))
#define THROW_ERROR_VALUE(text,value) throw runtime_error(global::to_string(boost::format(text) % value))
#define WRITE(text) cout << text << endl
#define WRITE_VALUE(text,value) cout << global::to_string(boost::format(text) % value) << endl
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
// global constants and functions used
//------------------------------------------------------------------------------
namespace global {
  static const int MPIROOT=0;

  template<typename T> string to_string(T value) {  //lexical_cast does funny things with double
    stringstream s1;

    s1 << value;
    return s1.str();
    }
//------------------------------------------------------------------------------
  template<typename T> T **make2DArray(int y,int x) {
    T **dest;

    dest=new T*[y];
    dest[0]=new T[y*x];
    for (int y1=1; y1<y; y1++)
      dest[y1]=&dest[0][y1*x];
    return dest;
    }
  //------------------------------------------------------------------------------
  string getFileName(string pathstring);
  }
//------------------------------------------------------------------------------
namespace CALC {
  #define IA 16807
  #define IM 2147483647
  #define AM (1.0/IM)
  #define IQ 127773
  #define IR 2836
  #define NTAB 32
  #define NDIV (1+(IM-1)/NTAB)
  #define EPS 1.2e-7
  #define RNMX (1.0-EPS)
  static long rseed=-123456789;

  void sran1(long rseed);
  double ran1();
  double Chi2(double x, int n);
  }
//------------------------------------------------------------------------------
#endif // GLOBAL_H
