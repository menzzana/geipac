#ifndef LOADER_H
#define LOADER_H
//------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <ctype.h>
#include "global.h"
//------------------------------------------------------------------------------
class Loader {
  public:
    static const char DELIMITER='\t';
    static const unsigned char ASCII0=48;

    string *splitDataString(string fstr, int ndatacolumns);
    int getColumnSize(string fstr);
    template<typename T> T *getEntry(T *first,string name);
    template<typename T> static T *addEntry(T *first);
    template<typename T> static T *loadFile(T *first,string filename);
  };
//------------------------------------------------------------------------------
class IMarkerData : public Loader {
  public:
    string markerid;
    IMarkerData *Next;

    IMarkerData();
    IMarkerData *getSingleRowData(string fstr);
    ~IMarkerData();
  };
//------------------------------------------------------------------------------
class IVariableData : public Loader {
  public:
    #define ENVIRONMENT "ENV"
    #define INDIVIDUAL_IDENTITY "INDID"
    string indid;
    double env,*cov;
    IVariableData *Next;

    IVariableData();
    IVariableData *getSingleRowData(string fstr);
    ~IVariableData();
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
    LimitData *getSingleRowData(string fstr);
    ~LimitData();
  };
//------------------------------------------------------------------------------
class FAMData : public Loader {
  public:
    enum genders {GENDER_UNKNOWN, MALE, FEMALE};
    enum phenotypes {PHENOTYPE_UNKNOWN, UNAFFECTED,AFFECTED};
    enum fam_file_position {FAM_FAMILYID,FAM_INDIVIDUALID,FAM_PATERNALID,FAM_MATERNAL_ID,FAM_GENDER,FAM_PHENOTYPE};
    static const int MAX_COLUMNS=6;
    #define FAM_FILE ".fam"

    int gender;
    int phenotype;
    string individualid;
    FAMData *Next;

    FAMData();
    FAMData *getSingleRowData(string fstr);
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
    BIMData *Next;

    BIMData();
    BIMData *getSingleRowData(string fstr);
    ~BIMData();
  };
//------------------------------------------------------------------------------
class BEDData : public Loader {
  public:
    #define BED_FILE ".bed"
    static const char HOMOZYGOTE1=0b00;
    static const char HOMOZYGOTE2=0b11;
    static const char HETEROZYGOTE=0b01;
    static const char MISSING_GENOTYPE=0b10;

    int genotype;
    FAMData *fam;
    BIMData *bim;
    BEDData *Next;

    BEDData();
    static BEDData *loadBinaryFile(string filename,FAMData *firstfam,BIMData *firstbim);
    ~BEDData();
  };
//------------------------------------------------------------------------------
#endif // LOADER_H
