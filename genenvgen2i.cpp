#include "genenvgen2i.h"
//------------------------------------------------------------------------------
using namespace GenEnvGen2I;
//---------------------------------------------------------------------------
Analysis::Analysis() {
  param.randomseed=CALC::rseed;
  param.appnegative=false;
  param.apcalculation=NO_PROPORTION;
  param.threshold=THRESHOLD;
  param.iterations=ITERATIONS;
  param.permutations=0;
  param.cutoff=CUTOFF;
  param.model=NO_MODEL;
  param.permutation_output=PERMUTATION_TOTALDATA;
  markerid=NULL;
  imarkerid=NULL;
  individualid=NULL;
  chromosome=NULL;
  imarkinteraction=NULL;
  interaction=NULL;
  interactionfromfile=NULL;
  covariate=NULL;
  cutoff_mult=NULL;
  cutoff_app=NULL;
  gender=NULL;
  phenotype=NULL;
  permphenotype=NULL;
  genotype=NULL;
  allele1=NULL;
  allele2=NULL;
  riskfactors=NULL;
  covdata1=NULL;
  covdata2=NULL;
  rephenotype=NULL;
  nimarkerid=0;
  nmarkerid=0;
  nindividualid=0;
  nlimit=0;
  ncovariate=0;
  }
//------------------------------------------------------------------------------
Analysis::~Analysis() {
  delete[] markerid;
  delete[] imarkerid;
  delete[] individualid;
  delete[] chromosome;
  if (interactionfromfile!=NULL)
    delete interactionfromfile;
  else
    delete imarkinteraction;
  delete interaction;
  delete[] covariate;
  delete cutoff_mult;
  delete cutoff_app;
  delete gender;
  delete phenotype;
  delete permphenotype;
  delete[] genotype;
  delete allele1;
  delete allele2;
  delete riskfactors;
  delete[] covdata1;
  delete[] covdata2;
  delete[] cleancovdata;
  delete rephenotype;
  delete cleanrephenotype;
  }
//------------------------------------------------------------------------------
void Analysis::initialize() {
  riskfactors=new int[nindividualid];
  interaction=new int[nindividualid];
  if (param.permutations>0)
    permphenotype=new int[nindividualid];
  covdata1=global::make2DArray<double>(nindividualid,MATRIX_INDEX_COV1+ncovariate);
  covdata2=global::make2DArray<double>(nindividualid,MATRIX_INDEX_COV2+ncovariate);
  cleancovdata=global::make2DArray<double>(nindividualid,MATRIX_INDEX_COV2+ncovariate);
  rephenotype=new double[nindividualid];
  cleanrephenotype=new double[nindividualid];
  for (int y1=0; y1<nindividualid; y1++) {
    rephenotype[y1]=phenotype[y1]-1;
    for (int x1=0; x1<ncovariate; x1++) {
      covdata1[y1][MATRIX_INDEX_COV1+x1]=covariate[y1][x1];
      covdata2[y1][MATRIX_INDEX_COV2+x1]=covariate[y1][x1];
      }
    }
  logreg1.createArrays(cleanrephenotype,cleancovdata,MATRIX_INDEX_COV2+ncovariate);
  logreg2.createArrays(cleanrephenotype,cleancovdata,MATRIX_INDEX_COV1+ncovariate);
  }
//------------------------------------------------------------------------------
void Analysis::run(int imarkeridx) {
  bool belowthreshold,riskabovecutoff;
  int recode;
  string results[RESULT_COLUMNS::LENGTH_RESULTS];
  char riskallele;

  CALC::sran1(param.randomseed);
  setInteraction(imarkeridx);
  for (int permidx=0; permidx<=param.permutations; permidx++) {
    if (permidx==0)
      WRITE(STATUS_TEXT::ORIGINAL_START);
    else
      WRITE_VALUE(STATUS_TEXT::PERMUTATION_START,permidx);
    for (int markeridx=0; markeridx<nmarkerid; markeridx++) {
      recode=0;
      memcpy(interaction,imarkinteraction,nindividualid*sizeof(int));
      for (int i1=0; i1<RESULT_COLUMNS::LENGTH_RESULTS; i1++)
        results[i1]="0.0";
      results[RESULT_COLUMNS::STABLELRM]="NA";
      results[RESULT_COLUMNS::STABLELRA]="NA";
      results[RESULT_COLUMNS::PERM]=global::to_string(param.permutations);
      results[RESULT_COLUMNS::INTERACTION]=getInteractionMarkerName(imarkeridx);
      results[RESULT_COLUMNS::CHR]=chromosome[markeridx];
      results[RESULT_COLUMNS::SNP]=markerid[markeridx];
      riskallele=calculateRiskAllele(markeridx,results);
      results[RESULT_COLUMNS::RISK]=riskallele;
      calculateRiskFactors(markeridx,riskallele,recode);
      riskabovecutoff=calculateRiskMatrix(NULL);
      logreg1.dimy=cleanData(markeridx,rephenotype,covdata2,MATRIX_INDEX_COV2+ncovariate);


      logreg1.clearArrays();
      belowthreshold=logreg1.gradientDescent(param.iterations,param.threshold);
      if (riskabovecutoff) {
        logreg2.dimy=cleanData(markeridx,rephenotype,covdata1,MATRIX_INDEX_COV1+ncovariate);
        logreg2.clearArrays();
        belowthreshold=logreg2.gradientDescent(param.iterations,param.threshold);
        logreg2.calculateZ();
        results[RESULT_COLUMNS::STABLELRM]=belowthreshold?STATUS_TEXT::POSITIVE_CONVERGENCE:STATUS_TEXT::NEGATIVE_CONVERGENCE;
        results[RESULT_COLUMNS::MULT]=global::to_string(logreg2.getMULTPropability(MATRIX_INDEX_A1mB1m));
        results[RESULT_COLUMNS::ORMIO]=global::to_string(logreg2.theta[MATRIX_INDEX_A1m]);
        results[RESULT_COLUMNS::ORMIOL]=global::to_string(logreg2.lowCI(MATRIX_INDEX_A1m));
        results[RESULT_COLUMNS::ORMIOH]=global::to_string(logreg2.highCI(MATRIX_INDEX_A1m));
        results[RESULT_COLUMNS::ORMIO]=global::to_string(logreg2.theta[MATRIX_INDEX_B1m]);
        results[RESULT_COLUMNS::ORMIOL]=global::to_string(logreg2.lowCI(MATRIX_INDEX_B1m));
        results[RESULT_COLUMNS::ORMIOH]=global::to_string(logreg2.highCI(MATRIX_INDEX_B1m));
        results[RESULT_COLUMNS::ORMIO]=global::to_string(logreg2.theta[MATRIX_INDEX_A1mB1m]);
        results[RESULT_COLUMNS::ORMIOL]=global::to_string(logreg2.lowCI(MATRIX_INDEX_A1mB1m));
        results[RESULT_COLUMNS::ORMIOH]=global::to_string(logreg2.highCI(MATRIX_INDEX_A1mB1m));


        }
      if (logreg1.theta[MATRIX_INDEX_A1B0]<0 &&
          logreg1.theta[MATRIX_INDEX_A1B0]<logreg1.theta[MATRIX_INDEX_A0B1] &&
          logreg1.theta[MATRIX_INDEX_A1B0]<logreg1.theta[MATRIX_INDEX_A1B1])
        recode=1;
      if (logreg1.theta[MATRIX_INDEX_A0B1]<0 &&
          logreg1.theta[MATRIX_INDEX_A0B1]<logreg1.theta[MATRIX_INDEX_A1B0] &&
          logreg1.theta[MATRIX_INDEX_A0B1]<logreg1.theta[MATRIX_INDEX_A1B1]) {
        recode=2;
        swapInteractions();
        }
      if (logreg1.theta[MATRIX_INDEX_A1B1]<0 &&
          logreg1.theta[MATRIX_INDEX_A1B1]<logreg1.theta[MATRIX_INDEX_A1B0] &&
          logreg1.theta[MATRIX_INDEX_A1B1]<logreg1.theta[MATRIX_INDEX_A0B1]) {
        recode=3;
        swapInteractions();
        }
      if (recode>0) {
        calculateRiskFactors(markeridx,riskallele,recode);
        riskabovecutoff=calculateRiskMatrix(results);
        }

      //printResults(cout,results);


      }



    shufflePhenotype();
    }

  }
//------------------------------------------------------------------------------
void Analysis::alleleSummaryCount(int *alleles,int markeridx) {
  int offset;

  for (int i1=0; i1<nindividualid; i1++) {
    offset=(phenotype[i1]==PHENOTYPE_AFFECTED?INDEX_CASE_PRIMARY:INDEX_CONTROL_PRIMARY);
    switch (genotype[i1][markeridx]) {
      case HOMOZYGOTE_PRIMARY:
        alleles[offset]+=2;
        break;
      case HOMOZYGOTE_SECONDARY:
        alleles[offset+INDEX_POSITION_OFFSET]+=2;
        break;
      case HETEROZYGOTE:
        alleles[offset]++;
        alleles[offset+INDEX_POSITION_OFFSET]++;
        break;
      }
    }
  }
//------------------------------------------------------------------------------
string Analysis::getInteractionMarkerName(int imarkeridx) {
  return markerid[imarkerid[imarkeridx]];
  }
//------------------------------------------------------------------------------
void Analysis::setInteraction(int imarkeridx) {
  double ratioriskalleleprimary,ratioriskallelesecondary;
  int riskhomozygote;
  int alleles[]={0,0,0,0};

  if (interactionfromfile!=NULL) {
    imarkinteraction=interactionfromfile;
    return;
    }
  if (imarkinteraction==NULL)
    imarkinteraction=new int[nindividualid];
  alleleSummaryCount(alleles,imarkeridx);
  ratioriskalleleprimary=(double)alleles[INDEX_CASE_PRIMARY]/(double)alleles[INDEX_CONTROL_PRIMARY];
  ratioriskallelesecondary=(double)alleles[INDEX_CASE_SECONDARY]/(double)alleles[INDEX_CONTROL_SECONDARY];
  riskhomozygote=(ratioriskalleleprimary>ratioriskallelesecondary?HOMOZYGOTE_PRIMARY:HOMOZYGOTE_SECONDARY);
  for (int i1=0; i1<nindividualid; i1++) {
    imarkinteraction[i1]=NA_INTERACTION;
    if (!validGeneticData(i1,imarkeridx))
      continue;
    switch (genotype[i1][imarkeridx]) {
      case HOMOZYGOTE_PRIMARY:
        imarkinteraction[i1]=(riskhomozygote==HOMOZYGOTE_PRIMARY?INTERACTION:NO_INTERACTION);
        break;
      case HOMOZYGOTE_SECONDARY:
        imarkinteraction[i1]=(riskhomozygote==HOMOZYGOTE_SECONDARY?INTERACTION:NO_INTERACTION);
        break;
      case HETEROZYGOTE:
        imarkinteraction[i1]=(isDominantOrXMale(i1,imarkeridx)?INTERACTION:NO_INTERACTION);
        break;
      }
    }
  }
//------------------------------------------------------------------------------
bool Analysis::validIndividualData(int individualidx,int markeridx) {
  return (genotype[individualidx][markeridx]!=ZYGOTE_UNKNOWN &&
          phenotype[individualidx]!=PHENOTYPE_UNKNOWN &&
          interaction[individualidx]!=NA_INTERACTION);
  }
//------------------------------------------------------------------------------
bool Analysis::validGeneticData(int individualidx,int markeridx) {
  return (genotype[individualidx][markeridx]!=ZYGOTE_UNKNOWN &&
          phenotype[individualidx]!=PHENOTYPE_UNKNOWN);
  }
//------------------------------------------------------------------------------
char Analysis::calculateRiskAllele(int markeridx, string *results) {
  char controlmaxallele,casemaxallele,caseminallele;
  int controlmax,casemax;
  int primarycount,secondarycount;
  double controlmaxratio,casemaxratio;
  int alleles[]={0,0,0,0};

  alleleSummaryCount(alleles,markeridx);
  primarycount=alleles[INDEX_CONTROL_PRIMARY]+alleles[INDEX_CASE_PRIMARY];
  secondarycount=alleles[INDEX_CONTROL_SECONDARY]+alleles[INDEX_CASE_SECONDARY];
  if (primarycount>secondarycount) {
    results[RESULT_COLUMNS::MAJOR]=allele1[markeridx];
    results[RESULT_COLUMNS::MINOR]=allele2[markeridx];
    }
  else {
    results[RESULT_COLUMNS::MAJOR]=allele2[markeridx];
    results[RESULT_COLUMNS::MINOR]=allele1[markeridx];
    }
  if (alleles[INDEX_CONTROL_PRIMARY]>alleles[INDEX_CONTROL_SECONDARY]) {
    controlmaxallele=allele1[markeridx];
    controlmax=alleles[INDEX_CONTROL_PRIMARY];
    }
  else {
    controlmaxallele=allele2[markeridx];
    controlmax=alleles[INDEX_CONTROL_SECONDARY];
    }
  if (alleles[INDEX_CASE_PRIMARY]>alleles[INDEX_CASE_SECONDARY]) {
    casemaxallele=allele1[markeridx];
    caseminallele=allele2[markeridx];
    casemax=alleles[INDEX_CASE_PRIMARY];
    }
  else {
    casemaxallele=allele2[markeridx];
    caseminallele=allele1[markeridx];
    casemax=alleles[INDEX_CASE_SECONDARY];
    }
  controlmaxratio=(double)controlmax/(double)(alleles[INDEX_CONTROL_PRIMARY]+alleles[INDEX_CONTROL_SECONDARY]);
  casemaxratio=(double)casemax/(double)(alleles[INDEX_CASE_PRIMARY]+alleles[INDEX_CASE_SECONDARY]);
  if (casemaxratio>controlmaxratio && casemaxallele==controlmaxallele)
    return casemaxallele;
  else
    return caseminallele;
  }
//------------------------------------------------------------------------------
bool Analysis::isDominantOrXMale(int individualidx,int markeridx) {
  if (param.model==DOMINANT)
    return true;
  return (boost::iequals(chromosome[markeridx],CHROMOSOME_X) && gender[individualidx]==GENDER_MALE);
  }
//------------------------------------------------------------------------------
void Analysis::calculateRiskFactors(int markeridx,char riskallele,int recode) {
  int healthygenotype,riskgenotype;

  healthygenotype=riskallele==allele1[markeridx]?HOMOZYGOTE_SECONDARY:HOMOZYGOTE_PRIMARY;
  riskgenotype=riskallele==allele1[markeridx]?HOMOZYGOTE_PRIMARY:HOMOZYGOTE_SECONDARY;
  for (int i1=0; i1<nindividualid; i1++) {
    riskfactors[i1]=NA_RISK;
    if (!validIndividualData(i1,markeridx))
      continue;
    if (genotype[i1][markeridx]==healthygenotype)
      riskfactors[i1]=NO_RISK;
    if (genotype[i1][markeridx]==riskgenotype)
      riskfactors[i1]=RISK;
    if (genotype[i1][markeridx]==HETEROZYGOTE)
      riskfactors[i1]=isDominantOrXMale(i1,markeridx)?RISK:NO_RISK;
    if (recode%2==1)
      riskfactors[i1]=(riskfactors[i1]==NO_RISK?RISK:NO_RISK);
    }
  }
//------------------------------------------------------------------------------
bool Analysis::calculateRiskMatrix(string *results) {
  static int NINDIVIDUALTYPES=8;
  int affstatus;

  int individualtypes[]={0,0,0,0,0,0,0,0};
  for (int y1=0; y1<nindividualid; y1++) {
    fill_n(covdata1[y1],MATRIX_INDEX_COV1,NA_INTERACTION);
    fill_n(covdata2[y1],MATRIX_INDEX_COV2,NA_INTERACTION);
    if (riskfactors[y1]==NA_RISK)
      continue;
    fill_n(covdata1[y1],MATRIX_INDEX_COV1,NO_INTERACTION);
    fill_n(covdata2[y1],MATRIX_INDEX_COV2,NO_INTERACTION);
    covdata1[y1][LRTHETA0]=1;
    covdata2[y1][LRTHETA0_M]=1;
    affstatus=rephenotype[y1];
    if (riskfactors[y1]==NO_RISK && interaction[y1]==NO_INTERACTION) {
      covdata2[y1][MATRIX_INDEX_A0B0]=INTERACTION;
      individualtypes[affstatus]++;
      }
    if (riskfactors[y1]==RISK && interaction[y1]==NO_INTERACTION) {
      covdata2[y1][MATRIX_INDEX_A1B0]=INTERACTION;
      individualtypes[RESULT_COLUMNS::IND10_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    if (riskfactors[y1]==NO_RISK && interaction[y1]>=INTERACTION) {
      covdata2[y1][MATRIX_INDEX_A0B1]=INTERACTION;
      individualtypes[RESULT_COLUMNS::IND01_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    if (riskfactors[y1]==RISK && interaction[y1]>=INTERACTION) {
      covdata2[y1][MATRIX_INDEX_A1B1]=INTERACTION;
      individualtypes[RESULT_COLUMNS::IND11_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    if (riskfactors[y1]==RISK)
      covdata1[y1][MATRIX_INDEX_A1m]=INTERACTION;
    if (interaction[y1]>=INTERACTION)
      covdata1[y1][MATRIX_INDEX_B1m]=INTERACTION;
    covdata1[y1][MATRIX_INDEX_A1mB1m]=covdata1[y1][MATRIX_INDEX_A1m]*covdata1[y1][MATRIX_INDEX_B1m];
    }
  if (results!=NULL) {
    results[RESULT_COLUMNS::IND00_0]=individualtypes[RESULT_COLUMNS::IND00_0-RESULT_COLUMNS::IND00_0];
    results[RESULT_COLUMNS::IND00_1]=individualtypes[RESULT_COLUMNS::IND00_1-RESULT_COLUMNS::IND00_0];
    results[RESULT_COLUMNS::IND01_0]=individualtypes[RESULT_COLUMNS::IND01_0-RESULT_COLUMNS::IND00_0];
    results[RESULT_COLUMNS::IND01_1]=individualtypes[RESULT_COLUMNS::IND01_1-RESULT_COLUMNS::IND00_0];
    results[RESULT_COLUMNS::IND10_0]=individualtypes[RESULT_COLUMNS::IND10_0-RESULT_COLUMNS::IND00_0];
    results[RESULT_COLUMNS::IND10_1]=individualtypes[RESULT_COLUMNS::IND10_1-RESULT_COLUMNS::IND00_0];
    results[RESULT_COLUMNS::IND11_0]=individualtypes[RESULT_COLUMNS::IND11_0-RESULT_COLUMNS::IND00_0];
    results[RESULT_COLUMNS::IND11_1]=individualtypes[RESULT_COLUMNS::IND11_1-RESULT_COLUMNS::IND00_0];
    }
  for (int i1=0; i1<NINDIVIDUALTYPES; i1++)
    if (individualtypes[i1]<=param.cutoff)
      return false;
  return true;
  }
//------------------------------------------------------------------------------
int Analysis::cleanData(int markeridx,double *y, double **x, int dimx) {
  int y2,x1;

  for (int y1=y2=0; y1<nindividualid; y1++) {
    if (!validIndividualData(y1,markeridx))
      continue;
    cleanrephenotype[y2]=y[y1];
    for (x1=0; x1<dimx; x1++) {
      if (x[y1][x1]==NA_INTERACTION)
        break;
      cleancovdata[y2][x1]=x[y1][x1];
      }
    if (x1==dimx)
      y2++;
    }
  return y2;
  }
//------------------------------------------------------------------------------
void Analysis::shufflePhenotype() {
  memcpy(permphenotype,phenotype,nindividualid*sizeof(int));
  for (int idx=0; idx<nindividualid; idx++)
    swap(permphenotype[idx],permphenotype[(int)CALC::ran1()*nindividualid]);
  }
//------------------------------------------------------------------------------
void Analysis::swapInteractions() {
  for (int i1=0; i1<nindividualid; i1++)
    if (interaction[i1]!=NA_INTERACTION)
      interaction[i1]=interaction[i1]==NO_INTERACTION?INTERACTION:NO_INTERACTION;
  }
//------------------------------------------------------------------------------
void Analysis::printResults(ostream &stream,string *results) {
  #define DELIMITER "\t "

  for (int i1=0; i1<RESULT_COLUMNS::LENGTH_RESULTS; i1++)
    if (results==NULL)
      stream << RESULT_COLUMNS::RESULT_COLUMN_TEXT[i1] << DELIMITER;
    else
      stream << results[i1] << DELIMITER;
  stream << endl;
  }
//------------------------------------------------------------------------------

