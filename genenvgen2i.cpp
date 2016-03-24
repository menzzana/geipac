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
  individualexist=NULL;
  chromosome=NULL;
  interaction=NULL;
  interactionfromfile=false;
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
  response1=NULL;
  response2=NULL;
  nimarkerid=0;
  nmarkerid=0;
  nindividualid=0;
  nlimit=0;
  ncovariate=0;
  }
//------------------------------------------------------------------------------
void Analysis::setInteractionFromData(int imarkeridx) {
  float ratioriskalleleprimary,ratioriskallelesecondary;
  int offset,riskhomozygote;
  float alleles[4]={0,0,0,0};

  if (interaction==NULL)
    interaction=new int[nindividualid];
  for (int i1=0; i1<nindividualid; i1++) {
    offset=(phenotype[i1]==PHENOTYPE_UNAFFECTED?INDEX_CASE_PRIMARY:INDEX_CONTROL_SECONDARY);
    switch (genotype[i1][imarkeridx]) {
      case HOMOZYGOTE_PRIMARY:
        alleles[offset]+=2;
        break;
      case HOMOZYGOTE_SECONDARY:
        alleles[offset+1]+=2;
        break;
      case HETEROZYGOTE:
        alleles[offset+1]++;
        alleles[offset]++;
        break;
      }
    }
  ratioriskalleleprimary=alleles[INDEX_CASE_PRIMARY]/alleles[INDEX_CONTROL_PRIMARY];
  ratioriskallelesecondary=alleles[INDEX_CASE_SECONDARY]/alleles[INDEX_CONTROL_SECONDARY];
  riskhomozygote=(ratioriskalleleprimary>ratioriskallelesecondary?HOMOZYGOTE_PRIMARY:HOMOZYGOTE_SECONDARY);
  for (int i1=0; i1<nindividualid; i1++) {
    switch (genotype[i1][imarkeridx]) {
      case HETEROZYGOTE:
        interaction[i1]=(param.model==DOMINANT?1:0);
        break;
      case HOMOZYGOTE_PRIMARY:
        interaction[i1]=(riskhomozygote==HOMOZYGOTE_PRIMARY?1:0);
        break;
      case HOMOZYGOTE_SECONDARY:
        interaction[i1]=(riskhomozygote==HOMOZYGOTE_SECONDARY?1:0);
        break;
      }
    }
	}
//------------------------------------------------------------------------------
void Analysis::initialize() {
  individualexist=new bool[nindividualid];
  riskfactors=new int[nindividualid];
  if (param.permutations>0)
    permphenotype=new int[nindividualid];
  if (ncovariate>0) {
    covdata1=global::make2DArray<double>(nindividualid,MATRIX_INDEX_COV1+ncovariate);
    covdata2=global::make2DArray<double>(nindividualid,MATRIX_INDEX_COV2+ncovariate);
    newcovdata=global::make2DArray<double>(nindividualid,MATRIX_INDEX_COV2+ncovariate);
    response1=new double[nindividualid];
    response2=new double[nindividualid];
    newresponse=new double[nindividualid];
    }
  }
//------------------------------------------------------------------------------
void Analysis::run(int imarkeridx) {
  LogReg *lr1,*lr2,*lr3;
  bool belowthreshold,riskabovecutoff;
  int recode,nnewdata;
  string results[RESULT_COLUMNS::LENGTH_RESULTS];
  char riskallele;

  CALC::sran1(param.randomseed);
  if (!interactionfromfile)
    setInteractionFromData(imarkeridx);
  for (int permidx=0; permidx<=param.permutations; permidx++) {
    if (permidx==0)
      WRITE(STATUS_TEXT::ORIGINAL_START);
    else
      WRITE_VALUE(STATUS_TEXT::PERMUTATION_START,permidx);
    for (int markeridx=0; markeridx<nmarkerid; markeridx++) {
      recode=0;
      for (int i1=0; i1<RESULT_COLUMNS::LENGTH_RESULTS; i1++)
        results[i1]="0.0";
      results[RESULT_COLUMNS::STABLELRM]="NA";
      results[RESULT_COLUMNS::STABLELRA]="NA";
      results[RESULT_COLUMNS::PERM]=param.permutations;
      results[RESULT_COLUMNS::INTERACTION]=imarkerid[imarkeridx];
      results[RESULT_COLUMNS::CHR]=chromosome[markeridx];
      results[RESULT_COLUMNS::SNP]=markerid[markeridx];
      setReducedDataMask(markeridx);
      riskallele=calculateRiskAllele(markeridx,results);
      results[RESULT_COLUMNS::RISK]=riskallele;
      calculateRiskFactors(markeridx,riskallele,recode);
      for (int x1=0; x1<ncovariate; x1++)
        for (int y1=0; y1<nindividualid; y1++) {
          covdata1[y1][MATRIX_INDEX_COV1+x1]=covariate[y1][x1];
          covdata2[y1][MATRIX_INDEX_COV2+x1]=covariate[y1][x1];
          }
      for (int i1=0; i1<nindividualid; i1++) {
        response1[i1]=NA_RISK;
        response2[i1]=NA_RISK;
        if (individualexist[i1]) {
          response1[i1]=phenotype[i1]-1;
          response2[i1]=phenotype[i1]-1;
          }
        }
      riskabovecutoff=calculateRiskMatrix(NULL);
      nnewdata=cleanData(response2,covdata2,MATRIX_INDEX_COV2+ncovariate);

      for (int i1=0; i1<nnewdata; i1++) {
        cout<<newresponse[i1]<<"-";

        for (int i2=0; i2<MATRIX_INDEX_COV2+ncovariate; i2++)
          cout<<newcovdata[i1][i2];
        cout<<endl;
        }


      lr1=new LogReg();
      lr1->createArrays(newresponse,newcovdata,nnewdata,ncovariate+MATRIX_INDEX_COV2);
      lr1->normalize();
      belowthreshold=lr1->gradientDescent(param.iterations,param.threshold);
      if (riskabovecutoff) {
        nnewdata=cleanData(response1,covdata1,MATRIX_INDEX_COV1+ncovariate);
        lr2=new LogReg();
        lr2->createArrays(newresponse,newcovdata,nnewdata,ncovariate+MATRIX_INDEX_COV1);
        lr2->clearArrays();
        lr2->normalize();
        belowthreshold=lr2->gradientDescent(param.iterations,param.threshold);
        lr2->calculateZ();
        results[RESULT_COLUMNS::STABLELRM]=belowthreshold?STATUS_TEXT::POSITIVE_CONVERGENCE:STATUS_TEXT::NEGATIVE_CONVERGENCE;
        results[RESULT_COLUMNS::MULT]=getMULTPropability(lr2,MATRIX_INDEX_A1mB1m+1);
        results[RESULT_COLUMNS::ORMIO]=lr2->theta[MATRIX_INDEX_A1m+1];
        results[RESULT_COLUMNS::ORMIOL]=lr2->lowCI(MATRIX_INDEX_A1m+1);
        results[RESULT_COLUMNS::ORMIOH]=lr2->highCI(MATRIX_INDEX_A1m+1);
        results[RESULT_COLUMNS::ORMIO]=lr2->theta[MATRIX_INDEX_B1m+1];
        results[RESULT_COLUMNS::ORMIOL]=lr2->lowCI(MATRIX_INDEX_B1m+1);
        results[RESULT_COLUMNS::ORMIOH]=lr2->highCI(MATRIX_INDEX_B1m+1);
        results[RESULT_COLUMNS::ORMIO]=lr2->theta[MATRIX_INDEX_A1mB1m+1];
        results[RESULT_COLUMNS::ORMIOL]=lr2->lowCI(MATRIX_INDEX_A1mB1m+1);
        results[RESULT_COLUMNS::ORMIOH]=lr2->highCI(MATRIX_INDEX_A1mB1m+1);

        }


      delete lr1;
      delete lr2;
      delete lr3;
      }



    shufflePhenotype();
    }

  }
//------------------------------------------------------------------------------
void Analysis::setReducedDataMask(int markeridx) {
  for (int i1=0; i1<nindividualid; i1++)
    individualexist[i1]=(genotype[i1][markeridx]!=ZYGOTE_UNKNOWN &&
                      phenotype[i1]!=PHENOTYPE_UNKNOWN &&
                      interaction[i1]!=MISSING_INTERACTION);
  }
//------------------------------------------------------------------------------
char Analysis::calculateRiskAllele(int markeridx, string *results) {
  char controlmaxallele,casemaxallele,caseminallele;
  int controlmax,casemax;
  int offset,primarycount,secondarycount;
  double controlmaxratio,casemaxratio;

  int alleles[]={0,0,0,0};
  for (int i1=0; i1<nindividualid; i1++)
    if (individualexist[i1]) {
      switch (phenotype[i1]) {
        case PHENOTYPE_UNAFFECTED:
          offset=INDEX_CONTROL_PRIMARY;
          break;
        case PHENOTYPE_AFFECTED:
          offset=INDEX_CASE_PRIMARY;
          break;
        }
      switch (genotype[i1][markeridx]) {
        case HOMOZYGOTE_PRIMARY:
          alleles[offset]+=2;
          break;
        case HETEROZYGOTE:
          alleles[offset]++;
          alleles[offset+INDEX_CONTROL_CASE_OFFSET]++;
          break;
        case HOMOZYGOTE_SECONDARY:
          alleles[offset+INDEX_CONTROL_CASE_OFFSET]+=2;
          break;
        }
      }
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
bool Analysis::isDominantOrXMale(int individualidx ) {
  if (param.model==DOMINANT)
    return true;
  return (boost::iequals(chromosome[individualidx],CHROMOSOME_X) && gender[individualidx]==GENDER_MALE);
  }
//------------------------------------------------------------------------------
void Analysis::calculateRiskFactors(int markeridx,char riskallele,int recode) {
  int healthygenotype,riskgenotype;

  healthygenotype=riskallele==allele1[markeridx]?HOMOZYGOTE_SECONDARY:HOMOZYGOTE_PRIMARY;
  riskgenotype=riskallele==allele1[markeridx]?HOMOZYGOTE_PRIMARY:HOMOZYGOTE_SECONDARY;
  for (int i1=0; i1<nindividualid; i1++) {
    riskfactors[i1]=NA_RISK;
    if (individualexist[i1]) {
      if (genotype[i1][markeridx]==healthygenotype)
        riskfactors[i1]=NO_RISK;
      if (genotype[i1][markeridx]==riskgenotype)
        riskfactors[i1]=RISK;
      if (genotype[i1][markeridx]==HETEROZYGOTE)
        riskfactors[i1]=isDominantOrXMale(i1)?RISK:NO_RISK;
      if (recode%2==1)
        riskfactors[i1]=(riskfactors[i1]==NO_RISK?RISK:NO_RISK);
      }
    }
  }
//------------------------------------------------------------------------------
bool Analysis::calculateRiskMatrix(string *results) {
  static int NINDIVIDUALTYPES=8;
  int affstatus;

  int individualtypes[]={0,0,0,0,0,0,0,0};
  for (int y1=0; y1<nindividualid; y1++)
    for (int x1=0; x1<MATRIX_INDEX_COV2; x1++) {
      if (x1<MATRIX_INDEX_COV1)
        covdata1[y1][x1]=NO_COVARIATE;
      covdata2[y1][x1]=NO_COVARIATE;
      }
  for (int y1=0; y1<nindividualid; y1++) {
    if (riskfactors[y1]==NA_RISK)
      continue;
    affstatus=phenotype[y1]-1;
    if (riskfactors[y1]==NO_RISK && interaction[y1]==NO_INTERACTION) {
      covdata2[y1][MATRIX_INDEX_A0B0]=INTERACTION;
      individualtypes[affstatus]++;
      }
    else
      covdata2[y1][MATRIX_INDEX_A0B0]=NO_INTERACTION;
    if (riskfactors[y1]==RISK && interaction[y1]==NO_INTERACTION) {
      covdata2[y1][MATRIX_INDEX_A1B0]=INTERACTION;
      individualtypes[RESULT_COLUMNS::IND10_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    else
      covdata2[y1][MATRIX_INDEX_A1B0]=NO_INTERACTION;
    if (riskfactors[y1]==NO_RISK && interaction[y1]>=INTERACTION) {
      covdata2[y1][MATRIX_INDEX_A0B1]=INTERACTION;
      individualtypes[RESULT_COLUMNS::IND01_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    else
      covdata2[y1][MATRIX_INDEX_A0B1]=NO_INTERACTION;
    if (riskfactors[y1]==RISK && interaction[y1]>=INTERACTION) {
      covdata2[y1][MATRIX_INDEX_A1B1]=INTERACTION;
      individualtypes[RESULT_COLUMNS::IND11_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    else
      covdata2[y1][MATRIX_INDEX_A0B1]=NO_INTERACTION;
    if (riskfactors[y1]==RISK)
      covdata1[y1][MATRIX_INDEX_A1m]=INTERACTION;
    else
      covdata1[y1][MATRIX_INDEX_A1m]=NO_INTERACTION;
    if (interaction[y1]>NO_INTERACTION)
      covdata1[y1][MATRIX_INDEX_B1m]=INTERACTION;
    else
      covdata1[y1][MATRIX_INDEX_A1m]=NO_INTERACTION;
    if (covdata1[y1][MATRIX_INDEX_A1m]!=MISSING_INTERACTION && covdata1[y1][MATRIX_INDEX_B1m]!=MISSING_INTERACTION)
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
int Analysis::cleanData(double *y, double **x, int dimx) {
  int y2,x1;

  for (int y1=y2=0; y1<nindividualid; y1++) {
    if (!individualexist[y1])
      continue;
    for (x1=0; x1<dimx; x1++)
      if (x[y1][x1]==NO_COVARIATE)
        break;
    if (x1==dimx) {
      for (int x1=0; x1<dimx; x1++)
        newcovdata[y2][x1]=x[y1][x1];
      newresponse[y2]=y[y1];
      y2++;
      }
    }
  return y2;
  }
//------------------------------------------------------------------------------
double Analysis::getMULTPropability(LogReg *lr,int idx) {
  return CALC::Chi2(lr->sumofsquares[idx]/pow(lr->z[idx],2),lr->dimy);
  }
//------------------------------------------------------------------------------
void Analysis::shufflePhenotype() {
  memcpy(permphenotype,phenotype,nindividualid*sizeof(int));
  for (int idx=0; idx<nindividualid; idx++)
    swap(permphenotype[idx],permphenotype[(int)CALC::ran1()*nindividualid]);
  }
//------------------------------------------------------------------------------
Analysis::~Analysis() {
  delete markerid;
  delete imarkerid;
  delete individualid;
  delete individualexist;
  delete chromosome;
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
  delete[] newcovdata;
  delete response1;
  delete response2;
  delete newresponse;
  }
//------------------------------------------------------------------------------

