#include "loader.h"
//------------------------------------------------------------------------------
template<typename T> T *Loader::getEntry(T *first,string name) {
  T *tl1;

  for (tl1=first; tl1!=NULL; tl1=tl1->Next)
    if (name.compare(tl1->markerid)==0)
      return tl1;
  return NULL;
  }
//------------------------------------------------------------------------------
template<typename T> T *Loader::addEntry(T *first) {
  T *tl1,*tl2;

  tl2=new T();
  if (first!=NULL) {
    for (tl1=first; tl1->Next!=NULL; tl1=tl1->Next);
    tl1->Next=tl2;
    }
  return tl2;
  }
//------------------------------------------------------------------------------
string *Loader::splitDataString(string fstr,int ndatacolumns) {
  int i1,i2;
  string *data;

  data=new string[ndatacolumns];
  for (i1=i2=0; i2<ndatacolumns; i2++) {
    data[i2]=fstr.substr(i1,fstr.find(DELIMITER,i1)-i1);
    i1=fstr.find(DELIMITER,i1)+1;
    if (i1==0)
      break;
    }
  return data;
  }
//---------------------------------------------------------------------------
int Loader::getColumnSize(string fstr) {
  int col,i1;

  for (col=i1=0; (i1=fstr.find(DELIMITER,i1)+1)>0; col++);
  return col+1;
  }
//---------------------------------------------------------------------------
template<typename T> T *Loader::loadFile(T *first,string filename) {
  ifstream fpr;
  string fstr;
  int nrows;
  T *data1,*data2;

  try {
    fpr.open(filename.c_str());
    if (!fpr.good())
      THROW_ERROR_VALUE(ERROR_TEXT::FILE_NOT_FOUND,filename);
    data1=data2=first;
    for (nrows=0; getline(fpr,fstr); nrows++) {
      data2=data2->getSingleRowData(fstr);
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
template IMarkerData *Loader::loadFile<IMarkerData>(IMarkerData *first,string filename);
template LimitData *Loader::loadFile<LimitData>(LimitData *first,string filename);
template IVariableData *Loader::loadFile<IVariableData>(IVariableData *first,string filename);
template FAMData *Loader::loadFile<FAMData>(FAMData *first,string filename);
template BIMData *Loader::loadFile<BIMData>(BIMData *first,string filename);
//---------------------------------------------------------------------------
IMarkerData::IMarkerData() {
  markerid="";
  Next=NULL;
  }
//---------------------------------------------------------------------------
IMarkerData *IMarkerData::getSingleRowData(string fstr) {
  IMarkerData *data1;

  data1=addEntry(this);
  data1->markerid=fstr;
  return data1;
  }
//---------------------------------------------------------------------
IMarkerData::~IMarkerData() {
  delete Next;
  }
//---------------------------------------------------------------------------
LimitData::LimitData() {
  cutoff_app=cutoff_mult=0;
  Next=NULL;
  }
//---------------------------------------------------------------------------
LimitData *LimitData::getSingleRowData(string fstr) {
  static int app_col=-1,mult_col=-1;
  LimitData *data1;
  string *splitdata;

  data1=NULL;
  splitdata=splitDataString(fstr,MAX_COLUMNS);
  if (app_col<0) {
    for (int i1=0; i1<MAX_COLUMNS; i1++) {
      if (boost::iequals(splitdata[0],CUTOFF_AP_PVALUE))
        app_col=i1;
      if (boost::iequals(splitdata[0],CUTOFF_MULTIPLICATIVE))
        mult_col=i1;
      }
    }
  else
    if (app_col>=0 && mult_col>=0) {
      data1=addEntry(this);
      data1->cutoff_app=atof(splitdata[app_col].c_str());
      data1->cutoff_mult=atof(splitdata[mult_col].c_str());
      }
  delete splitdata;
  return data1;
  }
//---------------------------------------------------------------------------
LimitData::~LimitData() {
  delete Next;
  }
//---------------------------------------------------------------------------
IVariableData::IVariableData() {
  cov=NULL;
  Next=NULL;
  }
//---------------------------------------------------------------------------
IVariableData *IVariableData::getSingleRowData(string fstr) {
  static int env_col=-1,ncol=0;
  static bool indid_col=false;
  IVariableData *data1;
  string *splitdata;

  data1=NULL;
  if (env_col<0)
    ncol=getColumnSize(fstr);
  splitdata=splitDataString(fstr,ncol);
  if (env_col<0) {
    indid_col=boost::iequals(splitdata[0],INDIVIDUAL_IDENTITY);
    for (int i1=1; i1<ncol; i1++)
      if (boost::iequals(splitdata[i1],ENVIRONMENT))
        env_col=i1;
    }
  else
    if (indid_col && env_col>=0) {
      data1=addEntry(this);
      data1->indid=splitdata[0];
      data1->env=atof(splitdata[env_col].c_str());
      data1->cov=new double[ncol-2];
      for (int i1=1; i1<ncol; i1++)
        data1->cov[i1<env_col?i1-1:i1-2]=atof(splitdata[i1].c_str());
      }
  delete splitdata;
  return data1;
  }
//---------------------------------------------------------------------------
IVariableData::~IVariableData() {
  delete cov;
  delete Next;
  }
//---------------------------------------------------------------------------
FAMData::FAMData() {
  gender=-1;
  phenotype=-1;
  individualid="";
  Next=NULL;
  }
//---------------------------------------------------------------------------
FAMData *FAMData::getSingleRowData(string fstr) {
  FAMData *data1;
  string *splitdata;

  splitdata=splitDataString(fstr,MAX_COLUMNS);
  data1=addEntry(this);
  data1->individualid=splitdata[FAM_INDIVIDUALID];
  data1->gender=atoi(splitdata[FAM_GENDER].c_str());
  data1->phenotype=atoi(splitdata[FAM_INDIVIDUALID].c_str());
  return data1;
  }
//---------------------------------------------------------------------------
FAMData::~FAMData() {
  delete Next;
  }
//---------------------------------------------------------------------------
BIMData::BIMData() {
  chromosome="";
  allele1=allele2=0;
  markerid="";
  Next=NULL;
  }
//---------------------------------------------------------------------------
BIMData *BIMData::getSingleRowData(string fstr) {
  BIMData *data1;
  string *splitdata;

  splitdata=splitDataString(fstr,MAX_COLUMNS);
  data1=addEntry(this);
  data1->chromosome=splitdata[BIM_CHROMOSOME];
  data1->markerid=splitdata[BIM_MARKER];
  data1->allele1=splitdata[BIM_ALLELE1].c_str()[0];
  data1->allele2=splitdata[BIM_ALLELE2].c_str()[0];
  return data1;
  }
//---------------------------------------------------------------------------
BIMData::~BIMData() {
  delete Next;
  }
//---------------------------------------------------------------------------
BEDData::BEDData() {
  genotype=-1;
  bim=NULL;
  fam=NULL;
  }
//---------------------------------------------------------------------------
BEDData *BEDData::loadBinaryFile(string filename,FAMData *firstfam,BIMData *firstbim) {
  static const char MAGIC1=0x6c;
  static const char MAGIC2=0x1b;

  char c1,c2;
  bool individual_major;
  ifstream fpr;
  BEDData *data1,*data2;
  BIMData *bim1;
  FAMData *fam1;

  try {
    fpr.open(filename.c_str());
    if (!fpr.good())
      return NULL;
    if (fpr.get()!=MAGIC1 && fpr.get()!=MAGIC2)
      THROW_ERROR(ERROR_TEXT::NO_BED_FILE);
    individual_major=fpr.get()>0;
    data1=NULL;
    fam1=firstfam;
    bim1=firstbim;
    while (!fpr.get(c1)) {
      for (int i1=0; i1<4; i1++) {
        if (bim1==NULL) {
          bim1=firstbim;
          if (!individual_major)
            break;
          }
        if (fam1==NULL) {
          fam1=firstfam;
          if (individual_major)
            break;
          }
        data2=addEntry(data2);
        if (data1==NULL)
          data1=data2;
        data2->fam=fam1;
        data2->bim=bim1;
        data2->genotype=(c1 & HOMOZYGOTE2);
        c1=(c1>>2);
        if (individual_major)
          bim1=bim1->Next;
        else
          fam1=fam1->Next;
        }
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
BEDData::~BEDData() {
  delete Next;
  }
//---------------------------------------------------------------------------
