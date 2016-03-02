#include "global.h"
//---------------------------------------------------------------------------
bool global::deleteResultFile(string filename) {
  if (!boost::filesystem::exists(filename))
    return false;
  boost::filesystem::remove(filename);
  return true;
  }
//---------------------------------------------------------------------------
string global::setOutputDirectory(string dirname) {
  string dirname1;
  int i1;

  dirname1=dirname;
  if (dirname1=="") {
    boost::posix_time::ptime now=boost::posix_time::second_clock::local_time();
    dirname1=FILE_TEXT::OUTPUT_DIRECTORY+to_string(now.date().day()+now.date().month()+now.date().year());
    for (i1=1; boost::filesystem::exists(dirname1+to_string(i1)); i1++);
    dirname1=dirname1+to_string(i1);
    }
  if (boost::filesystem::exists(dirname1))
    return dirname1;
  if (boost::filesystem::create_directory(dirname1))
    return dirname1;
  return "";
  }
//---------------------------------------------------------------------------
void CALC::sran1(long seedvalue) {
  // Initialize with negative number
  if (rseed<0)
    rseed=seedvalue;
  }
//---------------------------------------------------------------------------
double CALC::ran1() {
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (rseed<=0 || !iy) {
    if (-rseed<1)
      rseed=1;
    else
      rseed=-rseed;
    for (j=NTAB+7; j>=0; j--) {
      k=rseed/IQ;
      rseed=IA*(rseed-k*IQ)-IR*k;
      if (rseed<0)
        rseed+=IM;
      if (j<NTAB)
        iv[j]=rseed;
      }
    iy=iv[0];
    }
  k=(rseed)/IQ;
  rseed=IA*(rseed-k*IQ)-IR*k;
  if (rseed<0)
    rseed+=IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=rseed;
  temp=AM*iy;
  if (temp>RNMX)
    return RNMX;
  return temp;
  }
//---------------------------------------------------------------------------
