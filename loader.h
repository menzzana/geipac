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

#ifndef LOADER_H
#define LOADER_H
//------------------------------------------------------------------------------
#include <cstddef>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <ctype.h>
#include "global.h"
#include "genenvgen2i.h"
//------------------------------------------------------------------------------
class Loader {
  public:
    #define SEPARATOR "\t "
    static bool deleteResultFile(string filename);
    static string setOutputDirectory(string dirname);
    static string *splitDataString(string fstr, int ndatacolumns);
    static int getColumnSize(string fstr);
//------------------------------------------------------------------------------
// Class templates
//------------------------------------------------------------------------------
    template<typename T> T *getEntry(string name) {
      T *tl1;

      for (tl1=(T *)this; tl1!=nullptr; tl1=tl1->next)
        if (name.compare(tl1->markerid)==0)
          return tl1;
      return nullptr;
      }
//------------------------------------------------------------------------------
    template<typename T> T *push_back() {
      T *tl1,*tl2;

      try {
        tl2=new T();
        }
      catch(exception &e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
        }
      /*
       Check if this!=nullptr does not work in later GNU C++ (v6+)
       Therefore the -fno-delete-null-pointer-checks was added as compiler parameter
       See CMakeLists.txt
      */
      if (this!=nullptr) {
        for (tl1=(T *)this; tl1->next!=nullptr; tl1=tl1->next);
        tl1->next=tl2;
        }
      return tl2;
      }
//------------------------------------------------------------------------------
    template<typename T> int Length() {
      T *tl1;
      int i1;

      for (tl1=(T *)this,i1=0; tl1!=nullptr; tl1=tl1->next,i1++);
      return i1;
      }
//------------------------------------------------------------------------------
    template<typename T, typename K> T *get(T K::*pmember,int length, T *dest1) {
      K *tl1;
      T *dest;
      int i1;

      if (length==0)
        return nullptr;
      dest=dest1==nullptr?new T[length]:dest1;
      for (tl1=(K *)this,i1=0; tl1!=nullptr; tl1=tl1->next,i1++)
        dest[i1]=tl1->*pmember;
      return dest;
      }
//------------------------------------------------------------------------------
    template<typename T, typename K> T **get(T *K::*pmember,int rows,int columns, T **dest1) {
      K *tl1;
      T **dest;
      int y1,x1;

      if (rows==0 || columns==0)
        return nullptr;
      dest=dest1==nullptr?global::make2DArray<T>(rows,columns):dest1;
      for (tl1=(K *)this,y1=0; tl1!=nullptr; tl1=tl1->next,y1++)
        for (x1=0; x1<columns; x1++)
          dest[y1][x1]=(tl1->*pmember)[x1];
      return dest;
      }
//------------------------------------------------------------------------------
    template<typename T> static T *loadFile(string filename) {
      ifstream fpr;
      string fstr;
      T *data1,*data2;

      try {
        fpr.open(filename.c_str());
        if (!fpr.good())
          THROW_ERROR_VALUE(ERROR_TEXT::FILE_NOT_FOUND,filename);
        data1=data2=nullptr;
        while (getline(fpr,fstr)) {
          fstr=boost::algorithm::trim_copy(fstr);
          if (fstr.length()==0)
            continue;
          data2=data2->getSingleRowData(fstr,data1);
          if (data2==nullptr)
            continue;
          if (data1==nullptr)
            data1=data2;
          }
        fpr.close();
        if (data1==nullptr)
          THROW_ERROR_VALUE(ERROR_TEXT::NO_FILE_LOADED,filename);
        return data1;
        }
      catch(exception &e) {
        cerr << e.what() << endl;
        return nullptr;
        }
      }
//------------------------------------------------------------------------------
    template<typename T, typename K> bool areAllIndividualPresent(K *famdata) {
      K *fd1;
      T *dt1;

      if (this==nullptr)
        return true;
      for (dt1=(T *)this,fd1=famdata; dt1!=nullptr && fd1!=nullptr; dt1=dt1->next,fd1=fd1->next)
        if (dt1->individualid!=fd1->individualid) {
          WRITELN_VALUE(ERROR_TEXT::MISSING_INDIVIDUAL,dt1->individualid);
          return false;
          }
      return dt1==nullptr && fd1==nullptr;
      }
//------------------------------------------------------------------------------
    template<typename T> void Delete() {
      T *tl1,*tl2;

      tl1=(T *)this;
      while (tl1!=nullptr) {
        tl2=tl1;
        tl1=tl1->next;
        delete tl2;
        }
      }
  };
//------------------------------------------------------------------------------
class IMarkerData : public Loader {
  public:
    string markerid;
    int index;
    IMarkerData *next;

    IMarkerData();
    IMarkerData *getSingleRowData(string fstr,IMarkerData *first);
  };
//------------------------------------------------------------------------------
class LimitData : public Loader {
  public:
    #define CUTOFF_AP_PVALUE "CUTOFF_APP"
    #define CUTOFF_MULTIPLICATIVE "CUTOFF_MULT"
    static const int MAX_COLUMNS=2;
    double cutoff_mult,cutoff_app;
    LimitData *next;

    LimitData();
    LimitData *getSingleRowData(string fstr,...);
  };
//------------------------------------------------------------------------------
class FAMData : public Loader {
  public:
    enum POSITION {
      FAMILYID,INDIVIDUALID,PATERNALID,MATERNAL_ID,GENDER,PHENOTYPE,MAX_COLUMNS
      };
    #define FAM_FILE ".fam"

    int gender;
    int phenotype;
    string individualid;
    int index;
    FAMData *next;

    FAMData();
    FAMData *getSingleRowData(string fstr,...);
  };
//------------------------------------------------------------------------------
class BIMData : public Loader {
  public:
    enum POSITION {
      CHROMOSOME,MARKER,GENE_DISTANCE,BASE_POS,ALLELE1,ALLELE2,MAX_COLUMNS
      };
    #define BIM_FILE ".bim"

    char allele1,allele2;
    string markerid,chromosome;
    int index;
    BIMData *next;

    BIMData();
    BIMData *getSingleRowData(string fstr,...);
    bool setInteractionMarkerIndex(IMarkerData *imarker);
  };
//------------------------------------------------------------------------------
class BEDData : public Loader {
  public:
    #define BED_FILE ".bed"
    static const char ALL_BIT_SET=3;
    static const char MAGIC1=0x6c;
    static const char MAGIC2=0x1b;
    static const char GENOTYPES_PER_BYTE=4;
    static const char BITS_PER_GENOTYPE=2;

    int genotype;
    FAMData *fam;
    BIMData *bim;
    BEDData *next;

    BEDData();
    static BEDData *loadBinaryFile(string filename,FAMData *firstfam,BIMData *firstbim);
    int **getGenotypes(int y,int x);
  };
//------------------------------------------------------------------------------
class IVariableData : public Loader {
  public:
    #define ENVIRONMENT "ENV"
    #define INDIVIDUAL_IDENTITY "INDID"
    #define NA "NA"
    static const int ENV_NOVALUE=-1;
    string individualid;
    int interaction;
    IVariableData *next;

    IVariableData();
    IVariableData *getSingleRowData(string fstr,...);
    bool areInteractionsPresent();
  };
//------------------------------------------------------------------------------
class CovariateData : public Loader {
  public:
    #define INDIVIDUAL_IDENTITY "INDID"
    #define NA "NA"
    static const int COV_NOVALUE=0;
    string individualid;
    double *covariate;
    int ncovariate;
    CovariateData *next;

    CovariateData();
    ~CovariateData();
    CovariateData *getSingleRowData(string fstr,...);
  };
//------------------------------------------------------------------------------
class AltPhenotypeData : public Loader {
  public:
    enum POSITION {
      FAMILYID,INDIVIDUALID,PHENOTYPE1
      };
    string individualid;
    int *aphenotype,naphenotype;
    AltPhenotypeData *next;
    
    AltPhenotypeData();
    ~AltPhenotypeData();
    AltPhenotypeData *getSingleRowData(string fstr,...);  
  };
//------------------------------------------------------------------------------
#endif // LOADER_H
