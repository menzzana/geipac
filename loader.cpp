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
  fill_n(data,ndatacolumns,"");
  for (i1=i2=0; i2<ndatacolumns; i2++) {
    data[i2]=boost::algorithm::trim_copy(fstr.substr(i1,fstr.find_first_of(SEPARATOR,i1)-i1));
    i1=fstr.find_first_of(SEPARATOR,i1)+1;
    if (i1==0)
      break;
    }
  return data;
  }
//---------------------------------------------------------------------------
int Loader::getColumnSize(string fstr) {
  int col,i1;

  for (col=i1=0; (i1=fstr.find_first_of(SEPARATOR,i1)+1)>0; col++);
  return col+1;
  }
//---------------------------------------------------------------------------
IMarkerData::IMarkerData() {
  index=-1;
  markerid="";
  Next=nullptr;
  }
//---------------------------------------------------------------------------
IMarkerData *IMarkerData::getSingleRowData(string fstr,IMarkerData *first) {
  IMarkerData *data1;

  if (first->getEntry<IMarkerData>(fstr)!=nullptr)
    return nullptr;
  data1=this->push_back<IMarkerData>();
  data1->markerid=boost::algorithm::trim_copy(fstr);
  return data1;
  }
//---------------------------------------------------------------------------
LimitData::LimitData() {
  cutoff_app=cutoff_mult=0;
  Next=nullptr;
  }
//---------------------------------------------------------------------------
LimitData *LimitData::getSingleRowData(string fstr,...) {
  static int app_col=-1,mult_col=-1;
  LimitData *data1;
  string *splitdata;

  data1=nullptr;
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
      data1=this->push_back<LimitData>();
      data1->cutoff_app=atof(splitdata[app_col].c_str());
      data1->cutoff_mult=atof(splitdata[mult_col].c_str());
      }  
  delete[] splitdata;
  return data1;
  }
//---------------------------------------------------------------------------
FAMData::FAMData() {
  gender=(int)GenEnvGen2I::Gender::UNKNOWN;
  phenotype=(int)GenEnvGen2I::Phenotype::UNKNOWN;
  individualid="";
  index=0;
  Next=nullptr;
  }
//---------------------------------------------------------------------------
FAMData *FAMData::getSingleRowData(string fstr,...) {
  FAMData *data1;
  string *splitdata;
  static int index1=0;

  splitdata=splitDataString(fstr,MAX_COLUMNS);
  data1=this->push_back<FAMData>();
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
  Next=nullptr;
  }
//---------------------------------------------------------------------------
BIMData *BIMData::getSingleRowData(string fstr,...) {
  BIMData *data1;
  string *splitdata;
  static int index1=0;

  splitdata=splitDataString(fstr,MAX_COLUMNS);
  data1=this->push_back<BIMData>();
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

  if (imarker==nullptr)
    return true;
  for (imark1=imarker; imark1!=nullptr; imark1=imark1->Next) {
    for (bim1=this,idx=0; bim1!=nullptr; bim1=bim1->Next,idx++)
      if (boost::iequals(bim1->markerid,imark1->markerid))
        break;
    imark1->index=idx;
    if (bim1==nullptr) {
      WRITELN_VALUE(ERROR_TEXT::MISSING_MARKER,imark1->markerid);
      return false;
      }
    }
  return true;
  }
//---------------------------------------------------------------------------
BEDData::BEDData() {
  genotype=(int)GenEnvGen2I::Zygosity::UNKNOWN;
  bim=nullptr;
  fam=nullptr;
  Next=nullptr;
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
    data1=data2=nullptr;
    fam1=firstfam;
    bim1=firstbim;
    while (fpr.read(&c1,1)) {
      for (int i1=0; i1<GENOTYPES_PER_BYTE; i1++) {
        if (bim1==nullptr) {
          bim1=firstbim;
          fam1=fam1->Next;
          if (individual_major)
            break;
          }
        if (fam1==nullptr) {
          fam1=firstfam;
          bim1=bim1->Next;
          if (!individual_major)
            break;
          }
        data2=data2->push_back<BEDData>();
        if (data1==nullptr)
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
    return nullptr;
    }
  }
//---------------------------------------------------------------------------
int **BEDData::getGenotypes(int y,int x) {
  BEDData *bd1;
  int **genotypedest;

  if (y==0 || x==0)
    return nullptr;
  genotypedest=global::make2DArray<int>(y,x);
  for (bd1=this; bd1!=nullptr; bd1=bd1->Next)
    genotypedest[bd1->fam->index][bd1->bim->index]=bd1->genotype;
  return genotypedest;
  }
//---------------------------------------------------------------------------
IVariableData::IVariableData() {
  individualid="";
  interaction=ENV_NOVALUE;
  covariate=nullptr;
  ncovariate=0;
  Next=nullptr;
  }
//---------------------------------------------------------------------------
IVariableData::~IVariableData() {
  delete covariate;
  }
//---------------------------------------------------------------------------
IVariableData *IVariableData::getSingleRowData(string fstr,...) {
  static int env_col=-1,indid_col=-1,ncol=0;
  IVariableData *data1;
  string *splitdata;

  data1=nullptr;
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
    data1=this->push_back<IVariableData>();
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
bool IVariableData::areInteractionsPresent() {
  IVariableData *ivd1;

  for (ivd1=this; ivd1!=nullptr; ivd1=ivd1->Next)
    if (ivd1->interaction!=ENV_NOVALUE)
      return true;
  return false;
  }
//---------------------------------------------------------------------------
AltPhenotypeData::AltPhenotypeData() {
  naphenotype=0;
  aphenotype=nullptr;
  Next=nullptr;
  }
//---------------------------------------------------------------------------
AltPhenotypeData::~AltPhenotypeData() {
  delete aphenotype;
  }
//---------------------------------------------------------------------------
AltPhenotypeData *AltPhenotypeData::getSingleRowData(string fstr,...) {
  static int ncols=-1;
  AltPhenotypeData *data1;
  string *splitdata;

  data1=nullptr;
  if (ncols<0)
    ncols=getColumnSize(fstr);
  splitdata=splitDataString(fstr,ncols);
  data1=this->push_back<AltPhenotypeData>();
  data1->naphenotype=ncols-PHENOTYPE1;
  data1->individualid=splitdata[INDIVIDUALID];
  data1->aphenotype=new int[data1->naphenotype];
  for (int i1=PHENOTYPE1; i1<data1->naphenotype; i1++)
    data1->aphenotype[i1=i1-data1->naphenotype]=atoi(splitdata[i1].c_str());
  delete[] splitdata;
  return data1;
  }
//---------------------------------------------------------------------------
