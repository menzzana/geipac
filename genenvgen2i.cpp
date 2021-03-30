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

#include "genenvgen2i.h"
//------------------------------------------------------------------------------
using namespace GenEnvGen2I;
//------------------------------------------------------------------------------
Analysis::Analysis(DataStore *datastore) {
  data=datastore;
  imarkinteraction=nullptr;
  riskfactors=new RiskFactor[data->nindividualid];
  interaction=new int[data->nindividualid];
  covariate1=global::make2DArray<double>(data->nindividualid,(int)Apm::LENGTH+data->ncovariate);
  covariate2=global::make2DArray<double>(data->nindividualid,(int)Ap::LENGTH+data->ncovariate);
  for (int y1=0; y1<data->nindividualid; y1++)
    for (int x1=0; x1<data->ncovariate; x1++) {
      covariate1[y1][(int)Apm::LENGTH+x1]=data->covariate[y1][x1];
      covariate2[y1][(int)Ap::LENGTH+x1]=data->covariate[y1][x1];
      }  
  }
//------------------------------------------------------------------------------
Analysis::~Analysis() {
  if (data->interactionfromfile==nullptr)
    delete imarkinteraction;
  delete interaction;
  delete riskfactors;
  delete[] covariate1;
  delete[] covariate2;
  data=nullptr;
  }
//------------------------------------------------------------------------------
void Analysis::run(int interactivemarkeridx) {
  int npermuted[RESULT_COLUMNS::LENGTH_VALUES];
  string results_text[RESULT_COLUMNS::LENGTH_TEXT];
  double results_value[RESULT_COLUMNS::LENGTH_VALUES],original_value[RESULT_COLUMNS::LENGTH_VALUES];
  double perm_value[RESULT_COLUMNS::LENGTH_VALUES];

  setInteraction(interactivemarkeridx);
  #pragma omp for
    for (int markeridx=0; markeridx<data->nmarkerid; markeridx++) {
      for (int permidx=0; permidx<=data->permutations; permidx++) {
        if (permidx==0)
          WRITE(STATUS_TEXT::ORIGINAL_START);
        else
          WRITE(STATUS_TEXT::PERMUTATION_START);
        results_text[RESULT_COLUMNS::INTERACTION]=data->markerid[interactivemarkeridx];
        results_text[RESULT_COLUMNS::CHR]=data->chromosome[markeridx];
        results_text[RESULT_COLUMNS::SNP]=data->markerid[markeridx];
        results_text[RESULT_COLUMNS::PERM]=permidx==0?STATUS_TEXT::NO_PERMUTATION:global::to_string(permidx);
        analyzeData(markeridx,data->phenotype[permidx],results_text,results_value);
        if (permidx==0 || data->rawpermutation) {
          #pragma omp critical
            {
            printResults(*data->wres,results_text,RESULT_COLUMNS::LENGTH_TEXT);
            printResults(*data->wres,results_value,RESULT_COLUMNS::LENGTH_VALUES);
            *data->wres<<endl;
            }
          }
        if (data->permutations==0)
          continue;
        if (permidx==0) {
          fill_n(npermuted,RESULT_COLUMNS::LENGTH_VALUES,0);
          memcpy(original_value,results_value,RESULT_COLUMNS::LENGTH_VALUES*sizeof(double));
          perm_value[RESULT_COLUMNS::APP]=MAX_P_VALUE;
          perm_value[RESULT_COLUMNS::MULT]=MAX_P_VALUE;
          if (data->appnegative || results_value[RESULT_COLUMNS::APP]>0)
            perm_value[RESULT_COLUMNS::APP]=results_value[RESULT_COLUMNS::APP];
          if (results_value[RESULT_COLUMNS::MULT]>0)
            perm_value[RESULT_COLUMNS::MULT]=results_value[RESULT_COLUMNS::MULT];
          continue;
          }
        for (int residx=0; residx<RESULT_COLUMNS::LENGTH_VALUES; residx++)
          if (RESULT_COLUMNS::PERMUTED_VALUE[residx])
            switch(residx) {
              case RESULT_COLUMNS::STABLELRA:
              case RESULT_COLUMNS::STABLELRM:
                npermuted[residx]+=results_value[residx];
                break;
              case RESULT_COLUMNS::APP:
                if (data->appnegative || results_value[residx]>0) {
                  perm_value[residx]=min(perm_value[residx],results_value[residx]);
                  data->permuted_app[permidx]=min(data->permuted_app[permidx],results_value[residx]);
                  }
                break;
              case RESULT_COLUMNS::MULT:
                if (results_value[residx]>0) {
                  perm_value[residx]=min(perm_value[residx],results_value[residx]);
                  data->permuted_mult[permidx]=min(data->permuted_mult[permidx],results_value[residx]);
                  }
                break;
              default:
                if (results_value[residx]<=original_value[residx] && results_value[residx]>0)
                  npermuted[residx]++;
                break;
              }
        }
      if (data->permutations>0) {
        for (int residx=0; residx<RESULT_COLUMNS::LENGTH_VALUES; residx++)
          if (RESULT_COLUMNS::PERMUTED_VALUE[residx] && residx!=RESULT_COLUMNS::APP && residx!=RESULT_COLUMNS::MULT)
            perm_value[residx]=((double)npermuted[residx])/(double)data->permutations;
        #pragma omp critical
          {
          printResults(*data->wperm,results_text,RESULT_COLUMNS::PERM);
          printResults(*data->wperm,perm_value,RESULT_COLUMNS::LENGTH_VALUES,RESULT_COLUMNS::PERMUTED_VALUE);
          *data->wperm<<endl;
          }
        }
    }
  }
//------------------------------------------------------------------------------
void Analysis::analyzeData(int markeridx,int *phenotypex,string *results_text, double *results_value) {
  LogisticRegression logreg;
  bool belowthreshold;
  int recode,riskmatrix[N_RISK_MATRIX];
  char riskallele;

  boost::math::chi_squared chi2(1);
  recode=0;
  memcpy(interaction,imarkinteraction,data->nindividualid*sizeof(int));
  fill_n(results_value,RESULT_COLUMNS::LENGTH_VALUES,0);
  riskallele=calculateRiskAllele(markeridx,phenotypex,results_text);
  results_text[RESULT_COLUMNS::RISK]=riskallele;
  calculateRiskFactors(markeridx,phenotypex,riskallele,recode);
  calculateRiskMatrix(phenotypex,riskmatrix);
  if (!belowCutOff(riskmatrix)) {
    setCleanData(markeridx,phenotypex,covariate1,logreg.y,logreg.x,(int)Apm::LENGTH+data->ncovariate,true);
    belowthreshold=logreg.maximumLikelihoodRegression(data->iterations,data->threshold);
    results_value[RESULT_COLUMNS::STABLELRM]=belowthreshold?1:0;
    results_value[RESULT_COLUMNS::MULT]=1-cdf(chi2,pow(logreg.z((int)Lrm::A1B1),2));
    results_value[RESULT_COLUMNS::ORMIO]=logreg.oddsratio((int)Lrm::A1);
    results_value[RESULT_COLUMNS::ORMIOL]=logreg.lowCI((int)Lrm::A1);
    results_value[RESULT_COLUMNS::ORMIOH]=logreg.highCI((int)Lrm::A1);
    results_value[RESULT_COLUMNS::ORMOI]=logreg.oddsratio((int)Lrm::B1);
    results_value[RESULT_COLUMNS::ORMOIL]=logreg.lowCI((int)Lrm::B1);
    results_value[RESULT_COLUMNS::ORMOIH]=logreg.highCI((int)Lrm::B1);
    results_value[RESULT_COLUMNS::ORMII]=logreg.oddsratio((int)Lrm::A1B1);
    results_value[RESULT_COLUMNS::ORMIIL]=logreg.lowCI((int)Lrm::A1B1);
    results_value[RESULT_COLUMNS::ORMIIH]=logreg.highCI((int)Lrm::A1B1);
    double apmvalue;
    apmvalue=logreg.calculateAPMValue((int)Lrm::A1,(int)Lrm::B1,(int)Lrm::A1B1);
    results_value[RESULT_COLUMNS::APM]=apmvalue;
    double apmerror;
    apmerror=logreg.APSEM((int)Lrm::A1,(int)Lrm::B1,(int)Lrm::A1B1);
    if (apmerror>0) {
      results_value[RESULT_COLUMNS::APML]=logreg.lowCI(apmvalue,apmerror);
      results_value[RESULT_COLUMNS::APMH]=logreg.highCI(apmvalue,apmerror);
      boost::math::normal normaldist=boost::math::normal(0,apmerror);
      results_value[RESULT_COLUMNS::APMP]=(1-cdf(normaldist,abs(apmvalue)))*2;
      }
    }
  setCleanData(markeridx,phenotypex,covariate2,logreg.y,logreg.x,(int)Ap::LENGTH+data->ncovariate,false);
  belowthreshold=logreg.maximumLikelihoodRegression(data->iterations,data->threshold);
  if (logreg.beta((int)Lr::A1B0)<0 &&
      logreg.beta((int)Lr::A1B0)<logreg.beta((int)Lr::A0B1) &&
      logreg.beta((int)Lr::A1B0)<logreg.beta((int)Lr::A1B1))
    recode=1;
  if (logreg.beta((int)Lr::A0B1)<0 &&
      logreg.beta((int)Lr::A0B1)<logreg.beta((int)Lr::A1B0) &&
      logreg.beta((int)Lr::A0B1)<logreg.beta((int)Lr::A1B1)) {
    recode=2;
    swapInteractions();
    }
  if (logreg.beta((int)Lr::A1B1)<0 &&
      logreg.beta((int)Lr::A1B1)<logreg.beta((int)Lr::A1B0) &&
      logreg.beta((int)Lr::A1B1)<logreg.beta((int)Lr::A0B1)) {
    recode=3;
    swapInteractions();
    }
  if (recode>0) {
    calculateRiskFactors(markeridx,phenotypex,riskallele,recode);
    calculateRiskMatrix(phenotypex,riskmatrix);
    }
  results_value[RESULT_COLUMNS::IND00_0]=riskmatrix[RESULT_COLUMNS::IND00_0-RESULT_COLUMNS::IND00_0];
  results_value[RESULT_COLUMNS::IND00_1]=riskmatrix[RESULT_COLUMNS::IND00_1-RESULT_COLUMNS::IND00_0];
  results_value[RESULT_COLUMNS::IND01_0]=riskmatrix[RESULT_COLUMNS::IND01_0-RESULT_COLUMNS::IND00_0];
  results_value[RESULT_COLUMNS::IND01_1]=riskmatrix[RESULT_COLUMNS::IND01_1-RESULT_COLUMNS::IND00_0];
  results_value[RESULT_COLUMNS::IND10_0]=riskmatrix[RESULT_COLUMNS::IND10_0-RESULT_COLUMNS::IND00_0];
  results_value[RESULT_COLUMNS::IND10_1]=riskmatrix[RESULT_COLUMNS::IND10_1-RESULT_COLUMNS::IND00_0];
  results_value[RESULT_COLUMNS::IND11_0]=riskmatrix[RESULT_COLUMNS::IND11_0-RESULT_COLUMNS::IND00_0];
  results_value[RESULT_COLUMNS::IND11_1]=riskmatrix[RESULT_COLUMNS::IND11_1-RESULT_COLUMNS::IND00_0];
  results_value[RESULT_COLUMNS::RECODE]=recode;
  if (!belowCutOff(riskmatrix)) {
    setCleanData(markeridx,phenotypex,covariate2,logreg.y,logreg.x,(int)Ap::LENGTH+data->ncovariate,false);
    belowthreshold=logreg.maximumLikelihoodRegression(data->iterations,data->threshold);
    results_value[RESULT_COLUMNS::STABLELRA]=belowthreshold?1:0;
    results_value[RESULT_COLUMNS::ORIO]=logreg.oddsratio((int)Lr::A1B0);
    results_value[RESULT_COLUMNS::ORIOL]=logreg.lowCI((int)Lr::A1B0);
    results_value[RESULT_COLUMNS::ORIOH]=logreg.highCI((int)Lr::A1B0);
    results_value[RESULT_COLUMNS::ORII]=logreg.oddsratio((int)Lr::A1B1);
    results_value[RESULT_COLUMNS::ORIIL]=logreg.lowCI((int)Lr::A1B1);
    results_value[RESULT_COLUMNS::ORIIH]=logreg.highCI((int)Lr::A1B1);
    results_value[RESULT_COLUMNS::OROI]=logreg.oddsratio((int)Lr::A0B1);
    results_value[RESULT_COLUMNS::OROIL]=logreg.lowCI((int)Lr::A0B1);
    results_value[RESULT_COLUMNS::OROIH]=logreg.highCI((int)Lr::A0B1);
    double apvalue=0;
    switch(data->apcalculation) {
      case Proportion::DISEASE:
        apvalue=logreg.calculateRERI((int)Lr::A1B0,(int)Lr::A0B1,(int)Lr::A1B1)/logreg.oddsratio((int)Lr::A1B1);
        break;
      case Proportion::EFFECT:
        apvalue=logreg.calculateRERI((int)Lr::A1B0,(int)Lr::A0B1,(int)Lr::A1B1)/logreg.oddsratio((int)Lr::A1B1)-1;
        break;
      case Proportion::CORRECTED:
        apvalue=logreg.calculateRERI((int)Lr::A1B0,(int)Lr::A0B1,(int)Lr::A1B1)/
                max(logreg.oddsratio((int)Lr::A1B1),logreg.oddsratio((int)Lr::A1B0)+logreg.oddsratio((int)Lr::A0B1)-1);
        break;
      }
    results_value[RESULT_COLUMNS::AP]=apvalue;
    double aperror;
    aperror=logreg.APSEM((int)Lr::A1B0,(int)Lr::A0B1,(int)Lr::A1B1);
    if (aperror>0) {
      results_value[RESULT_COLUMNS::APL]=logreg.lowCI(apvalue,aperror);
      results_value[RESULT_COLUMNS::APH]=logreg.highCI(apvalue,aperror);
      boost::math::normal normaldist=boost::math::normal(0,aperror);
      results_value[RESULT_COLUMNS::APP]=(1-cdf(normaldist,abs(apvalue)))*2;
      }
    }  
  }
//------------------------------------------------------------------------------
void Analysis::printTotalPermutation(DataStore &data1) {
  int napp,nmult;

  if (!data1.totalPermutations() || data1.permutations==0)
    return;
  for (int sigidx=0; sigidx<data1.nlimit; sigidx++) {
    napp=nmult=0;
    for (int permidx=1; permidx<=data1.permutations; permidx++) {
      if (data1.permuted_app[permidx]<=data1.cutoff_app[sigidx])
        napp++;
      if (data1.permuted_mult[permidx]<=data1.cutoff_mult[sigidx])
        nmult++;
      }
    if (data1.cutoff_app[sigidx]>0)
      *data1.wtotperm<<DELIMITER<<data1.cutoff_app[sigidx]<<DELIMITER<<(double)napp/(double)data1.permutations;
    else
      *data1.wtotperm<<DELIMITER<<DELIMITER;
    if (data1.cutoff_mult[sigidx]>0)
      *data1.wtotperm<<DELIMITER<<data1.cutoff_mult[sigidx]<<DELIMITER<<(double)nmult/(double)data1.permutations;
    *data1.wtotperm<<endl;
    }
  for (int permidx=1; permidx<=data1.permutations; permidx++)
    *data1.wtotperm<<permidx<<DELIMITER<<DELIMITER<<data1.permuted_app[permidx]<<DELIMITER<<DELIMITER<<data1.permuted_mult[permidx]<<endl;
  }
//------------------------------------------------------------------------------
void Analysis::alleleSummaryCount(int *alleles,int markeridx,int *phenotypex) {
  int offset;

  for (int i1=0; i1<data->nindividualid; i1++) {
    offset=(phenotypex[i1]==(int)Phenotype::AFFECTED?(int)Allele::CASE_PRIMARY:(int)Allele::CONTROL_PRIMARY);
    switch ((Zygosity)data->genotype[i1][markeridx]) {
      case Zygosity::HOMOZYGOTE_PRIMARY:
        alleles[offset]+=2;
        break;
      case Zygosity::HOMOZYGOTE_SECONDARY:
        alleles[offset+(int)Allele::OFFSET]+=2;
        break;
      case Zygosity::HETEROZYGOTE:
        alleles[offset]++;
        alleles[offset+(int)Allele::OFFSET]++;
        break;
      }
    }
  }
//------------------------------------------------------------------------------
void Analysis::setInteraction(int interactivemarkeridx) {
  double ratioriskalleleprimary,ratioriskallelesecondary;
  int riskhomozygote;
  int alleles[]={0,0,0,0};

  if (data->interactionfromfile!=nullptr) {
    imarkinteraction=data->interactionfromfile;
    return;
    }
  if (imarkinteraction==nullptr)
    imarkinteraction=new int[data->nindividualid];
  alleleSummaryCount(alleles,interactivemarkeridx,data->phenotype[ORIGINAL]);
  ratioriskalleleprimary=(double)alleles[(int)Allele::CASE_PRIMARY]/(double)alleles[(int)Allele::CONTROL_PRIMARY];
  ratioriskallelesecondary=(double)alleles[(int)Allele::CASE_SECONDARY]/(double)alleles[(int)Allele::CONTROL_SECONDARY];
  riskhomozygote=(ratioriskalleleprimary>ratioriskallelesecondary?(int)Zygosity::HOMOZYGOTE_PRIMARY:(int)Zygosity::HOMOZYGOTE_SECONDARY);
  for (int i1=0; i1<data->nindividualid; i1++) {
    imarkinteraction[i1]=(int)Interaction::NA;
    if (!validGeneticData(i1,interactivemarkeridx,data->phenotype[ORIGINAL]))
      continue;
    switch ((Zygosity)data->genotype[i1][interactivemarkeridx]) {
      case Zygosity::HOMOZYGOTE_PRIMARY:
        imarkinteraction[i1]=(int)(riskhomozygote==(int)Zygosity::HOMOZYGOTE_PRIMARY?Interaction::YES:Interaction::NO);
        break;
      case Zygosity::HOMOZYGOTE_SECONDARY:
        imarkinteraction[i1]=(int)(riskhomozygote==(int)Zygosity::HOMOZYGOTE_SECONDARY?Interaction::YES:Interaction::NO);
        break;
      case Zygosity::HETEROZYGOTE:
        imarkinteraction[i1]=(int)(isDominantOrXMale(i1,interactivemarkeridx)?Interaction::YES:Interaction::NO);
        break;
      }
    }
  }
//------------------------------------------------------------------------------
bool Analysis::validIndividualData(int individualidx,int markeridx,int *phenotypex) {
  return (data->genotype[individualidx][markeridx]!=(int)Zygosity::UNKNOWN &&
          phenotypex[individualidx]!=(int)Phenotype::UNKNOWN &&
          interaction[individualidx]!=(int)Interaction::NA);
  }
//------------------------------------------------------------------------------
bool Analysis::validGeneticData(int individualidx,int markeridx,int *phenotypex) {
  return (data->genotype[individualidx][markeridx]!=(int)Zygosity::UNKNOWN &&
          phenotypex[individualidx]!=(int)Phenotype::UNKNOWN);
  }
//------------------------------------------------------------------------------
char Analysis::calculateRiskAllele(int markeridx, int *phenotypex, string *results) {
  char controlmaxallele,casemaxallele,caseminallele;
  int controlmax,casemax;
  int primarycount,secondarycount;
  double controlmaxratio,casemaxratio;
  int alleles[]={0,0,0,0};

  alleleSummaryCount(alleles,markeridx,phenotypex);
  primarycount=alleles[(int)Allele::CONTROL_PRIMARY]+alleles[(int)Allele::CASE_PRIMARY];
  secondarycount=alleles[(int)Allele::CONTROL_SECONDARY]+alleles[(int)Allele::CASE_SECONDARY];
  if (primarycount>secondarycount) {
    results[RESULT_COLUMNS::MAJOR]=data->allele1[markeridx];
    results[RESULT_COLUMNS::MINOR]=data->allele2[markeridx];
    }
  else {
    results[RESULT_COLUMNS::MAJOR]=data->allele2[markeridx];
    results[RESULT_COLUMNS::MINOR]=data->allele1[markeridx];
    }
  if (alleles[(int)Allele::CONTROL_PRIMARY]>alleles[(int)Allele::CONTROL_SECONDARY]) {
    controlmaxallele=data->allele1[markeridx];
    controlmax=alleles[(int)Allele::CONTROL_PRIMARY];
    }
  else {
    controlmaxallele=data->allele2[markeridx];
    controlmax=alleles[(int)Allele::CONTROL_SECONDARY];
    }
  if (alleles[(int)Allele::CASE_PRIMARY]>alleles[(int)Allele::CASE_SECONDARY]) {
    casemaxallele=data->allele1[markeridx];
    caseminallele=data->allele2[markeridx];
    casemax=alleles[(int)Allele::CASE_PRIMARY];
    }
  else {
    casemaxallele=data->allele2[markeridx];
    caseminallele=data->allele1[markeridx];
    casemax=alleles[(int)Allele::CASE_SECONDARY];
    }
  controlmaxratio=(double)controlmax/(double)(alleles[(int)Allele::CONTROL_PRIMARY]+alleles[(int)Allele::CONTROL_SECONDARY]);
  casemaxratio=(double)casemax/(double)(alleles[(int)Allele::CASE_PRIMARY]+alleles[(int)Allele::CASE_SECONDARY]);
  if (casemaxratio>controlmaxratio && casemaxallele==controlmaxallele)
    return casemaxallele;
  else
    return caseminallele;
  }
//------------------------------------------------------------------------------
bool Analysis::isDominantOrXMale(int individualidx,int markeridx) {
  if (data->model==Model::DOMINANT)
    return true;
  return (boost::iequals(data->chromosome[markeridx],CHROMOSOME_X) && data->gender[individualidx]==(int)Gender::MALE);
  }
//------------------------------------------------------------------------------
void Analysis::calculateRiskFactors(int markeridx,int *phenotypex,char riskallele,int recode) {
  int healthygenotype,riskgenotype;

  healthygenotype=riskallele==data->allele1[markeridx]?(int)Zygosity::HOMOZYGOTE_SECONDARY:(int)Zygosity::HOMOZYGOTE_PRIMARY;
  riskgenotype=riskallele==data->allele1[markeridx]?(int)Zygosity::HOMOZYGOTE_PRIMARY:(int)Zygosity::HOMOZYGOTE_SECONDARY;
  for (int i1=0; i1<data->nindividualid; i1++) {
    riskfactors[i1]=RiskFactor::NA;
    if (!validIndividualData(i1,markeridx,phenotypex))
      continue;
    if (data->genotype[i1][markeridx]==healthygenotype)
      riskfactors[i1]=RiskFactor::NO;
    if (data->genotype[i1][markeridx]==riskgenotype)
      riskfactors[i1]=RiskFactor::YES;
    if (data->genotype[i1][markeridx]==(int)Zygosity::HETEROZYGOTE)
      riskfactors[i1]=isDominantOrXMale(i1,markeridx)?RiskFactor::YES:RiskFactor::NO;
    if (recode%2==1)
      riskfactors[i1]=riskfactors[i1]==RiskFactor::NO?RiskFactor::YES:RiskFactor::NO;
    }
  }
//------------------------------------------------------------------------------
void Analysis::calculateRiskMatrix(int *phenotypex,int *riskmatrix) {
  int affstatus;

  fill_n(riskmatrix,N_RISK_MATRIX,0);
  for (int y1=0; y1<data->nindividualid; y1++) {
    fill_n(covariate1[y1],(int)Apm::LENGTH,(int)Interaction::NA);
    fill_n(covariate2[y1],(int)Ap::LENGTH,(int)Interaction::NA);
    if (riskfactors[y1]==RiskFactor::NA)
      continue;
    fill_n(covariate1[y1],(int)Apm::LENGTH,(int)Interaction::NO);
    fill_n(covariate2[y1],(int)Ap::LENGTH,(int)Interaction::NO);
    affstatus=phenotypex[y1]-1;
    if (riskfactors[y1]==RiskFactor::NO && interaction[y1]==(int)Interaction::NO) {
      covariate2[y1][(int)Ap::A0B0]=(int)Interaction::YES;
      riskmatrix[affstatus]++;
      }
    if (riskfactors[y1]==RiskFactor::YES && interaction[y1]==(int)Interaction::NO) {
      covariate2[y1][(int)Ap::A1B0]=(int)Interaction::YES;
      riskmatrix[RESULT_COLUMNS::IND10_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    if (riskfactors[y1]==RiskFactor::NO && interaction[y1]>=(int)Interaction::YES) {
      covariate2[y1][(int)Ap::A0B1]=(int)Interaction::YES;
      riskmatrix[RESULT_COLUMNS::IND01_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    if (riskfactors[y1]==RiskFactor::YES && interaction[y1]>=(int)Interaction::YES) {
      covariate2[y1][(int)Ap::A1B1]=(int)Interaction::YES;
      riskmatrix[RESULT_COLUMNS::IND11_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    if (riskfactors[y1]==RiskFactor::YES)
      covariate1[y1][(int)Apm::A1]=(int)Interaction::YES;
    if (interaction[y1]>=(int)Interaction::YES)
      covariate1[y1][(int)Apm::B1]=(int)Interaction::YES;
    covariate1[y1][(int)Apm::A1B1]=covariate1[y1][(int)Apm::A1]*covariate1[y1][(int)Apm::B1];
    }
  }
//------------------------------------------------------------------------------
bool Analysis::belowCutOff(int *riskmatrix) {
  for (int i1=0; i1<N_RISK_MATRIX; i1++)
    if (riskmatrix[i1]<=data->cutoff)
      return true;
  return false;
  }
//------------------------------------------------------------------------------
void Analysis::setCleanData(int markeridx,int *y, double **x, VectorXd &desty, MatrixXd &destx, int dimx, bool a0b0) {
  int y2,x1,x2;

  desty.resize(data->nindividualid);
  destx.resize(data->nindividualid,dimx+1-!a0b0);
  for (int y1=y2=0; y1<data->nindividualid; y1++) {
    if (!validIndividualData(y1,markeridx,y))
      continue;
    desty(y2)=y[y1]-1;
    destx(y2,(int)Lr::BETA0)=1;
    for (x1=0,x2=1; x1<dimx; x1++) {
      if (x[y1][x1]==(int)Interaction::NA)
        break;
      if (x1==(int)Ap::A0B0 && !a0b0)
        continue;
      destx(y2,x2)=x[y1][x1];
      x2++;
      }
    if (x1==dimx)
      y2++;
    }
  desty.conservativeResize(y2);
  destx.conservativeResize(y2,NoChange_t());
  }
//------------------------------------------------------------------------------
void Analysis::swapInteractions() {
  for (int i1=0; i1<data->nindividualid; i1++)
    if (interaction[i1]!=(int)Interaction::NA)
      interaction[i1]=(int)(interaction[i1]==(int)Interaction::NO?Interaction::YES:Interaction::NO);
  }
//------------------------------------------------------------------------------


