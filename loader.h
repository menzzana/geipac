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
    #define DELIMITER "\t "

    static bool deleteResultFile(string filename);
    static string setOutputDirectory(string dirname);
    string *splitDataString(string fstr, int ndatacolumns);
    int getColumnSize(string fstr);
//------------------------------------------------------------------------------
// Class templates
//------------------------------------------------------------------------------
    template<typename T> T *getEntry(string name) {
      T *tl1;

      for (tl1=(T *)this; tl1!=NULL; tl1=tl1->Next)
        if (name.compare(tl1->markerid)==0)
          return tl1;
      return NULL;
      }
//------------------------------------------------------------------------------
    template<typename T> T *addEntry() {
      T *tl1,*tl2;

      tl2=new T();
      if (this!=NULL) {
        for (tl1=(T *)this; tl1->Next!=NULL; tl1=tl1->Next);
        tl1->Next=tl2;
        }
      return tl2;
      }
//------------------------------------------------------------------------------
    template<typename T> int getLength() {
      T *tl1;
      int i1;

      for (tl1=(T *)this,i1=1; tl1->Next!=NULL; tl1=tl1->Next,i1++);
      return i1;
      }
//------------------------------------------------------------------------------
    template<typename T, typename K> T* get(T K::*pmember,int length) {
      K *tl1;
      T *dest;
      int i1;

      if (length==0)
        return NULL;
      dest=new T[length];
      for (tl1=(K *)this,i1=0; tl1!=NULL; tl1=tl1->Next,i1++)
        dest[i1]=tl1->*pmember;
      return dest;
      }
//---------------------------------------------------------------------------
    template<typename T> T *loadFile(string filename) {
      ifstream fpr;
      string fstr;
      int nrows;
      T *data1,*data2;

      try {
        fpr.open(filename.c_str());
        if (!fpr.good())
          THROW_ERROR_VALUE(ERROR_TEXT::FILE_NOT_FOUND,filename);
        data1=data2=(T *)this;
        for (nrows=0; getline(fpr,fstr); nrows++) {
          data2=data2->getSingleRowData(fstr,data1);
          if (data2==NULL)
            continue;
          if (data1==NULL)
            data1=data2;
          }
        fpr.close();
        return data1;
        }
      catch(exception &e) {
        cerr << e.what() << endl;
        return NULL;
        }
      }
//---------------------------------------------------------------------------
  };
//------------------------------------------------------------------------------
class IMarkerData : public Loader {
  public:
    string markerid;
    IMarkerData *Next;

    IMarkerData();
    IMarkerData *getSingleRowData(string fstr,IMarkerData *first);
    ~IMarkerData();
  };
//------------------------------------------------------------------------------
class LimitData : public Loader {
  public:
    #define CUTOFF_AP_PVALUE "CUTOFF_APP"
    #define CUTOFF_MULTIPLICATIVE "CUTOFF_MULT"
    static const int MAX_COLUMNS=2;
    double cutoff_mult,cutoff_app;
    LimitData *Next;

    LimitData();
    LimitData *getSingleRowData(string fstr,...);
    ~LimitData();
  };
//------------------------------------------------------------------------------
class FAMData : public Loader {
  public:
    enum fam_file_position {FAM_FAMILYID,FAM_INDIVIDUALID,FAM_PATERNALID,FAM_MATERNAL_ID,FAM_GENDER,FAM_PHENOTYPE};
    static const int MAX_COLUMNS=6;
    #define FAM_FILE ".fam"

    int gender;
    int phenotype;
    string individualid;
    int index;
    FAMData *Next;

    FAMData();
    FAMData *getSingleRowData(string fstr,...);
    ~FAMData();
  };
//------------------------------------------------------------------------------
class BIMData : public Loader {
  public:
    enum bim_file_position {BIM_CHROMOSOME,BIM_MARKER,BIM_GENE_DISTANCE,BIM_BASE_POS,BIM_ALLELE1,BIM_ALLELE2};
    static const int MAX_COLUMNS=6;
    #define BIM_FILE ".bim"

    char allele1,allele2;
    string markerid,chromosome;
    int index;
    BIMData *Next;

    BIMData();
    BIMData *getSingleRowData(string fstr,...);
    bool areInteractionMarkersPresent(IMarkerData *imarker);
    ~BIMData();
  };
//------------------------------------------------------------------------------
class BEDData : public Loader {
  public:
    #define BED_FILE ".bed"
    static const char ALL_BIT_SET=3;

    int genotype;
    FAMData *fam;
    BIMData *bim;
    BEDData *Next;

    BEDData();
    static BEDData *loadBinaryFile(string filename,FAMData *firstfam,BIMData *firstbim);
    int **getGenotypes(int y,int x);
    ~BEDData();
  };
//------------------------------------------------------------------------------
class IVariableData : public Loader {
  public:
    #define ENVIRONMENT "ENV"
    #define INDIVIDUAL_IDENTITY "INDID"
    #define NA "NA"
    static const int ENV_NOVALUE=-1;
    static const int COV_NOVALUE=0;
    string individualid;
    int interaction,*covariate;
    IVariableData *Next;

    IVariableData();
    IVariableData *getSingleRowData(string fstr,...);
    bool areAllIndividualPresent(FAMData *famdata);
    int **getCovariates(int y,int x);
    bool areInteractionsPresent();
    ~IVariableData();
  };
//------------------------------------------------------------------------------
#endif // LOADER_H
