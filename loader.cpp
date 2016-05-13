#include "loader.h"
//------------------------------------------------------------------------------
bool Loader::deleteResultFile(string filename) {
  if (!boost::filesystem::exists(filename))
    return false;
  boost::filesystem::remove(filename);
  return true;
  }
//---------------------------------------------------------------------------
string Loader::setOutputDirectory(string dirname) {
  string dirname1;
  int i1;

  dirname1=dirname;
  if (dirname1=="") {
    boost::posix_time::ptime now=boost::posix_time::second_clock::local_time();
    dirname1=FILE_TEXT::OUTPUT_DIRECTORY+global::to_string(now.date().day()+now.date().month()+now.date().year());
    for (i1=1; boost::filesystem::exists(dirname1+global::to_string(i1)); i1++);
    dirname1=dirname1+global::to_string(i1);
    }
  if (boost::filesystem::exists(dirname1))
    return dirname1;
  if (boost::filesystem::create_directory(dirname1))
    return dirname1;
  return "";
  }
//---------------------------------------------------------------------------
string *Loader::splitDataString(string fstr,int ndatacolumns) {
  int i1,i2;
  string *data;

  data=new string[ndatacolumns];
  for (i1=i2=0; i2<ndatacolumns; i2++) {
    data[i2]=boost::algorithm::trim_copy(fstr.substr(i1,fstr.find_first_of(DELIMITER,i1)-i1));
    i1=fstr.find_first_of(DELIMITER,i1)+1;
    if (i1==0)
      break;
    }
  return data;
  }
//---------------------------------------------------------------------------
int Loader::getColumnSize(string fstr) {
  int col,i1;

  for (col=i1=0; (i1=fstr.find_first_of(DELIMITER,i1)+1)>0; col++);
  return col+1;
  }
//---------------------------------------------------------------------------
IMarkerData::IMarkerData() {
  index=-1;
  markerid="";
  Next=NULL;
  }
//---------------------------------------------------------------------------
IMarkerData::~IMarkerData() {
  delete Next;
  }
//---------------------------------------------------------------------------
IMarkerData *IMarkerData::getSingleRowData(string fstr,IMarkerData *first) {
  IMarkerData *data1;

  if (first->getEntry<IMarkerData>(fstr)!=NULL)
    return NULL;
  data1=this->addEntry<IMarkerData>();
  data1->markerid=boost::algorithm::trim_copy(fstr);
  return data1;
  }
//---------------------------------------------------------------------------
LimitData::LimitData() {
  cutoff_app=cutoff_mult=0;
  Next=NULL;
  }
//---------------------------------------------------------------------------
LimitData::~LimitData() {
  delete Next;
  }
//---------------------------------------------------------------------------
LimitData *LimitData::getSingleRowData(string fstr,...) {
  static int app_col=-1,mult_col=-1;
  LimitData *data1;
  string *splitdata;

  data1=NULL;
  splitdata=splitDataString(fstr,MAX_COLUMNS);
  if (app_col<0) {
    for (int i1=0; i1<MAX_COLUMNS; i1++) {
      if (boost::iequals(splitdata[i1],CUTOFF_AP_PVALUE))
        app_col=i1;
      if (boost::iequals(splitdata[i1],CUTOFF_MULTIPLICATIVE))
        mult_col=i1;
      }
    }
  else
    if (app_col>=0 && mult_col>=0) {
      data1=this->addEntry<LimitData>();
      data1->cutoff_app=atof(splitdata[app_col].c_str());
      data1->cutoff_mult=atof(splitdata[mult_col].c_str());
      }
  delete[] splitdata;
  return data1;
  }
//---------------------------------------------------------------------------
FAMData::FAMData() {
  gender=GenEnvGen2I::GENDER_UNKNOWN;
  phenotype=GenEnvGen2I::PHENOTYPE_UNKNOWN;
  individualid="";
  index=0;
  Next=NULL;
  }
//---------------------------------------------------------------------------
FAMData::~FAMData() {
  delete Next;
  }
//---------------------------------------------------------------------------
FAMData *FAMData::getSingleRowData(string fstr,...) {
  FAMData *data1;
  string *splitdata;
  static int index1=0;

  splitdata=splitDataString(fstr,MAX_COLUMNS);
  data1=this->addEntry<FAMData>();
  data1->individualid=splitdata[POSITION::INDIVIDUALID];
  data1->gender=atoi(splitdata[POSITION::GENDER].c_str());
  data1->phenotype=atoi(splitdata[POSITION::PHENOTYPE].c_str());
  data1->index=index1;
  index1++;
  delete[] splitdata;
  return data1;
  }
//---------------------------------------------------------------------------
BIMData::BIMData() {
  chromosome="";
  allele1=allele2=index=0;
  markerid="";
  Next=NULL;
  }
//---------------------------------------------------------------------------
BIMData::~BIMData() {
  delete Next;
  }
//---------------------------------------------------------------------------
BIMData *BIMData::getSingleRowData(string fstr,...) {
  BIMData *data1;
  string *splitdata;
  static int index1=0;

  splitdata=splitDataString(fstr,MAX_COLUMNS);
  data1=this->addEntry<BIMData>();
  data1->chromosome=splitdata[POSITION::CHROMOSOME];
  data1->markerid=splitdata[POSITION::MARKER];
  data1->allele1=splitdata[POSITION::ALLELE1].c_str()[0];
  data1->allele2=splitdata[POSITION::ALLELE2].c_str()[0];
  data1->index=index1;
  index1++;
  delete[] splitdata;
  return data1;
  }
//---------------------------------------------------------------------------
bool BIMData::setInteractionMarkerIndex(IMarkerData *imarker) {
  IMarkerData *imark1;
  BIMData *bim1;
  int idx;

  if (imarker==NULL)
    return true;
  for (imark1=imarker; imark1!=NULL; imark1=imark1->Next) {
    for (bim1=this,idx=0; bim1!=NULL; bim1=bim1->Next,idx++)
      if (boost::iequals(bim1->markerid,imark1->markerid))
        break;
    imark1->index=idx;
    if (bim1==NULL) {
      WRITE_VALUE(ERROR_TEXT::MISSING_MARKER,imark1->markerid);
      return false;
      }
    }
  return true;
  }
//---------------------------------------------------------------------------
BEDData::BEDData() {
  genotype=GenEnvGen2I::ZYGOTE_UNKNOWN;
  bim=NULL;
  fam=NULL;
  Next=NULL;
  }
//---------------------------------------------------------------------------
BEDData::~BEDData() {
  delete Next;
  }
//---------------------------------------------------------------------------
BEDData *BEDData::loadBinaryFile(string filename,FAMData *firstfam,BIMData *firstbim) {
  char c1,c2;
  bool individual_major;
  ifstream fpr;
  BEDData *data1,*data2;
  BIMData *bim1;
  FAMData *fam1;

  try {
    fpr.open(filename.c_str(),ios::binary);
    if (!fpr.good())
      THROW_ERROR_VALUE(ERROR_TEXT::FILE_NOT_FOUND,filename);
    fpr.read(&c1,1);
    fpr.read(&c2,1);
    if (c1!=MAGIC1 && c2!=MAGIC2)
      THROW_ERROR(ERROR_TEXT::NO_BED_FILE);
    fpr.read(&c1,1);
    individual_major=(c1==0);
    data1=data2=NULL;
    fam1=firstfam;
    bim1=firstbim;
    while (fpr.read(&c1,1)) {
      for (int i1=0; i1<GENOTYPES_PER_BYTE; i1++) {
        if (bim1==NULL) {
          bim1=firstbim;
          fam1=fam1->Next;
          if (individual_major)
            break;
          }
        if (fam1==NULL) {
          fam1=firstfam;
          bim1=bim1->Next;
          if (!individual_major)
            break;
          }
        data2=data2->addEntry<BEDData>();
        if (data1==NULL)
          data1=data2;
        data2->fam=fam1;
        data2->bim=bim1;
        data2->genotype=(c1 & ALL_BIT_SET);
        c1=(c1>>BITS_PER_GENOTYPE);
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
int **BEDData::getGenotypes(int y,int x) {
  BEDData *bd1;
  int **genotypedest;

  if (y==0 || x==0)
    return NULL;
  genotypedest=global::make2DArray<int>(y,x);
  for (bd1=this; bd1!=NULL; bd1=bd1->Next)
    genotypedest[bd1->fam->index][bd1->bim->index]=bd1->genotype;
  return genotypedest;
  }
//---------------------------------------------------------------------------
IVariableData::IVariableData() {
  individualid="";
  interaction=ENV_NOVALUE;
  covariate=NULL;
  ncovariate=0;
  Next=NULL;
  }
//---------------------------------------------------------------------------
IVariableData::~IVariableData() {
  delete covariate;
  delete Next;
  }
//---------------------------------------------------------------------------
IVariableData *IVariableData::getSingleRowData(string fstr,...) {
  static int env_col=-1,indid_col=-1,ncol=0;
  IVariableData *data1;
  string *splitdata;

  data1=NULL;
  if (indid_col<0)
    ncol=getColumnSize(fstr);
  splitdata=splitDataString(fstr,ncol);
  if (indid_col<0) {
    for (int i1=1; i1<ncol; i1++) {
      if (boost::iequals(splitdata[i1],ENVIRONMENT))
        env_col=i1;
      if (boost::iequals(splitdata[i1],INDIVIDUAL_IDENTITY))
        indid_col=i1;
      }
    }
  else {
    data1=this->addEntry<IVariableData>();
    data1->individualid=splitdata[indid_col];
    if (env_col>=0)
      if (!boost::iequals(splitdata[env_col].c_str(),NA))
        data1->interaction=atoi(splitdata[env_col].c_str());
    data1->ncovariate=ncol-(env_col<0?1:2);
    data1->covariate=new int[data1->ncovariate];
    for (int i1=0,i2=0; i1<ncol; i1++) {
      if (i1==env_col || i1==indid_col)
        continue;
      if (boost::iequals(splitdata[i1].c_str(),NA))
        data1->covariate[i2]=COV_NOVALUE;
      else
        data1->covariate[i2]=atoi(splitdata[i1].c_str());
      i2++;
      }
    }
  delete[] splitdata;
  return data1;
  }
//---------------------------------------------------------------------------
bool IVariableData::areAllIndividualPresent(FAMData *famdata) {
  FAMData *fd1;
  IVariableData *ivd1;

  for (ivd1=this,fd1=famdata; ivd1!=NULL && fd1!=NULL; ivd1=ivd1->Next,fd1=fd1->Next)
    if (ivd1->individualid!=fd1->individualid) {
      WRITE_VALUE(ERROR_TEXT::MISSING_INDIVIDUAL,ivd1->individualid);
      return false;
      }
  return ivd1==NULL && fd1==NULL;
  }
//---------------------------------------------------------------------------
bool IVariableData::areInteractionsPresent() {
  IVariableData *ivd1;

  for (ivd1=this; ivd1!=NULL; ivd1=ivd1->Next)
    if (ivd1->interaction!=ENV_NOVALUE)
      return true;
  return false;
  }
//---------------------------------------------------------------------------
