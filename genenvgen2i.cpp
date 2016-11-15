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
  param.rawpermutation=false;
  param.totalpermutation=false;
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
  permuted_mult=NULL;
  permuted_app=NULL;
  gender=NULL;
  phenotype=NULL;
  genotype=NULL;
  allele1=NULL;
  allele2=NULL;
  riskfactors=NULL;
  covariate1=NULL;
  covariate2=NULL;
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
  delete permuted_mult;
  delete permuted_app;
  delete gender;
  delete[] phenotype;
  delete[] genotype;
  delete allele1;
  delete allele2;
  delete riskfactors;
  delete[] covariate1;
  delete[] covariate2;
  }
//------------------------------------------------------------------------------
void Analysis::initialize() {
  riskfactors=new int[nindividualid];
  interaction=new int[nindividualid];
  phenotype=global::make2DArray<int>(param.permutations+1,nindividualid);
  if (param.permutations>0) {
    permuted_mult=new double[param.permutations+1];
    permuted_app=new double[param.permutations+1];
    fill_n(permuted_mult,param.permutations+1,MAX_P_VALUE);
    fill_n(permuted_app,param.permutations+1,MAX_P_VALUE);
    }
  }
//------------------------------------------------------------------------------
void Analysis::createCovariateMatrix() {
  covariate1=global::make2DArray<double>(nindividualid,MATRIX_INDEX_COV1+ncovariate);
  covariate2=global::make2DArray<double>(nindividualid,MATRIX_INDEX_COV2+ncovariate);
  for (int y1=0; y1<nindividualid; y1++)
    for (int x1=0; x1<ncovariate; x1++) {
      covariate1[y1][MATRIX_INDEX_COV1+x1]=covariate[y1][x1];
      covariate2[y1][MATRIX_INDEX_COV2+x1]=covariate[y1][x1];
      }
  }
//------------------------------------------------------------------------------
void Analysis::permutePhenotypes() {
  CALC::sran1(param.randomseed);
  for (int permidx=1; permidx<=param.permutations; permidx++) {
    memcpy(phenotype[permidx],phenotype[permidx-1],nindividualid*sizeof(int));
    CALC::randomShuffle(phenotype[permidx],nindividualid);
    }
  }
//------------------------------------------------------------------------------
void Analysis::run(int interactivemarkeridx) {
  int npermuted[RESULT_COLUMNS::LENGTH_VALUES];
  string results_text[RESULT_COLUMNS::LENGTH_TEXT];
  double results_value[RESULT_COLUMNS::LENGTH_VALUES],original_value[RESULT_COLUMNS::LENGTH_VALUES];
  double perm_value[RESULT_COLUMNS::LENGTH_VALUES];
  
  setInteraction(interactivemarkeridx);
  for (int markeridx=0; markeridx<nmarkerid; markeridx++) {
    for (int permidx=0; permidx<=param.permutations; permidx++) {
      if (permidx==0)
        WRITE(STATUS_TEXT::ORIGINAL_START);
      else
        WRITE(STATUS_TEXT::PERMUTATION_START);
      results_text[RESULT_COLUMNS::INTERACTION]=markerid[interactivemarkeridx];
      results_text[RESULT_COLUMNS::CHR]=chromosome[markeridx];
      results_text[RESULT_COLUMNS::SNP]=markerid[markeridx];
      results_text[RESULT_COLUMNS::PERM]=permidx==0?STATUS_TEXT::NO_PERMUTATION:global::to_string(permidx);
      analyzeData(markeridx,phenotype[permidx],results_text,results_value);      
      if (permidx==0 || param.rawpermutation) {
        printResults(*param.wres,results_text,RESULT_COLUMNS::LENGTH_TEXT);
        printResults(*param.wres,results_value,RESULT_COLUMNS::LENGTH_VALUES);
        *param.wres<<endl;
        }
      if (param.permutations==0)
        continue;
      if (permidx==0) {
        fill_n(npermuted,RESULT_COLUMNS::LENGTH_VALUES,0);
        memcpy(original_value,results_value,RESULT_COLUMNS::LENGTH_VALUES*sizeof(double));
        perm_value[RESULT_COLUMNS::APP]=MAX_P_VALUE;
        perm_value[RESULT_COLUMNS::MULT]=MAX_P_VALUE;
        if (param.appnegative || results_value[RESULT_COLUMNS::APP]>0)
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
              if (param.appnegative || results_value[residx]>0) {
                perm_value[residx]=min(perm_value[residx],results_value[residx]);
                permuted_app[permidx]=min(permuted_app[permidx],results_value[residx]);
                }
              break;
            case RESULT_COLUMNS::MULT:
              if (results_value[residx]>0) {
                perm_value[residx]=min(perm_value[residx],results_value[residx]);
                permuted_mult[permidx]=min(permuted_mult[permidx],results_value[residx]);
                }
              break;
            default:
              if (results_value[residx]<=original_value[residx] && results_value[residx]>0)
                npermuted[residx]++;
              break;
            }
      }
    if (param.permutations>0) {
      printResults(*param.wperm,results_text,RESULT_COLUMNS::PERM);
      for (int residx=0; residx<RESULT_COLUMNS::LENGTH_VALUES; residx++)
        if (RESULT_COLUMNS::PERMUTED_VALUE[residx] && residx!=RESULT_COLUMNS::APP && residx!=RESULT_COLUMNS::MULT)
          perm_value[residx]=((double)npermuted[residx])/(double)param.permutations;
      printResults(*param.wperm,perm_value,RESULT_COLUMNS::LENGTH_VALUES,RESULT_COLUMNS::PERMUTED_VALUE);
      *param.wperm<<endl;
      }
    }  
  if (param.totalpermutation && param.permutations>0) {
    int napp,nmult;
    *param.wtotperm<<TOTAL_PERMUTATION;
    for (int sigidx=0; sigidx<nlimit; sigidx++) {
      napp=nmult=0;
      for (int permidx=1; permidx<=param.permutations; permidx++) {
        if (permuted_app[permidx]<=cutoff_app[sigidx])
          napp++;
        if (permuted_mult[permidx]<=cutoff_mult[sigidx])
          nmult++;
        }
      if (cutoff_app[sigidx]>0)
        *param.wtotperm<<DELIMITER<<cutoff_app[sigidx]<<DELIMITER<<(double)napp/(double)param.permutations;
      else
        *param.wtotperm<<DELIMITER<<DELIMITER;
      if (cutoff_mult[sigidx]>0)
        *param.wtotperm<<DELIMITER<<cutoff_mult[sigidx]<<DELIMITER<<(double)nmult/(double)param.permutations;
      *param.wtotperm<<endl;
      }
    for (int permidx=1; permidx<=param.permutations; permidx++)
      *param.wtotperm<<permidx<<DELIMITER<<DELIMITER<<permuted_app[permidx]<<DELIMITER<<DELIMITER<<permuted_mult[permidx]<<endl;
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
  memcpy(interaction,imarkinteraction,nindividualid*sizeof(int));
  fill_n(results_value,RESULT_COLUMNS::LENGTH_VALUES,0);
  riskallele=calculateRiskAllele(markeridx,phenotypex,results_text);
  results_text[RESULT_COLUMNS::RISK]=riskallele;
  calculateRiskFactors(markeridx,phenotypex,riskallele,recode);
  calculateRiskMatrix(phenotypex,riskmatrix);
  setCleanData(markeridx,phenotypex,covariate2,logreg.y,logreg.x,MATRIX_INDEX_COV2+ncovariate,false);
  belowthreshold=logreg.maximumLikelihoodRegression(param.iterations,param.threshold);

  cout<<"LR ->"<<markerid[markeridx]<<" ->";
  for (int i1=0; i1<logreg.beta.rows(); i1++)
    cout<<logreg.beta(i1)<<"; ";
  cout<<endl;

  if (!belowCutOff(riskmatrix)) {
    setCleanData(markeridx,phenotypex,covariate1,logreg.y,logreg.x,MATRIX_INDEX_COV1+ncovariate,true);
    belowthreshold=logreg.maximumLikelihoodRegression(param.iterations,param.threshold);

    cout<<"LR1 ->"<<markerid[markeridx]<<" ->";
    for (int i1=0; i1<logreg.beta.rows(); i1++)
      cout<<logreg.beta(i1)<<"; ";
    cout<<endl;


    results_value[RESULT_COLUMNS::STABLELRM]=belowthreshold?1:0;
    results_value[RESULT_COLUMNS::MULT]=cdf(chi2,pow(logreg.z(LR_INDEX_A1mB1m),2));
    results_value[RESULT_COLUMNS::ORMIO]=logreg.oddsratio(LR_INDEX_A1m);
    results_value[RESULT_COLUMNS::ORMIOL]=logreg.lowCI(LR_INDEX_A1m);
    results_value[RESULT_COLUMNS::ORMIOH]=logreg.highCI(LR_INDEX_A1m);
    results_value[RESULT_COLUMNS::ORMOI]=logreg.oddsratio(LR_INDEX_B1m);
    results_value[RESULT_COLUMNS::ORMOIL]=logreg.lowCI(LR_INDEX_B1m);
    results_value[RESULT_COLUMNS::ORMOIH]=logreg.highCI(LR_INDEX_B1m);
    results_value[RESULT_COLUMNS::ORMII]=logreg.oddsratio(LR_INDEX_A1mB1m);
    results_value[RESULT_COLUMNS::ORMIIL]=logreg.lowCI(LR_INDEX_A1mB1m);
    results_value[RESULT_COLUMNS::ORMIIH]=logreg.highCI(LR_INDEX_A1mB1m);
    double apmvalue;
    apmvalue=logreg.calculateAPMValue(LR_INDEX_A1m,LR_INDEX_B1m,LR_INDEX_A1mB1m);
    results_value[RESULT_COLUMNS::APM]=apmvalue;
    double apmerror;
    apmerror=logreg.APSEM(LR_INDEX_A1m,LR_INDEX_B1m,LR_INDEX_A1mB1m);
    results_value[RESULT_COLUMNS::APML]=logreg.lowCI(apmvalue,apmerror);
    results_value[RESULT_COLUMNS::APMH]=logreg.highCI(apmvalue,apmerror);
    boost::math::normal normaldist=boost::math::normal(0,apmerror);
    results_value[RESULT_COLUMNS::APMP]=1-cdf(normaldist,abs(apmvalue))*2;
    }
  if (logreg.beta(LR_INDEX_A1B0)<0 &&
      logreg.beta(LR_INDEX_A1B0)<logreg.beta(LR_INDEX_A0B1) &&
      logreg.beta(LR_INDEX_A1B0)<logreg.beta(LR_INDEX_A1B1))
    recode=1;
  if (logreg.beta(LR_INDEX_A0B1)<0 &&
      logreg.beta(LR_INDEX_A0B1)<logreg.beta(LR_INDEX_A1B0) &&
      logreg.beta(LR_INDEX_A0B1)<logreg.beta(LR_INDEX_A1B1)) {
    recode=2;
    swapInteractions();
    }
  if (logreg.beta(LR_INDEX_A1B1)<0 &&
      logreg.beta(LR_INDEX_A1B1)<logreg.beta(LR_INDEX_A1B0) &&
      logreg.beta(LR_INDEX_A1B1)<logreg.beta(LR_INDEX_A0B1)) {
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
    setCleanData(markeridx,phenotypex,covariate2,logreg.y,logreg.x,MATRIX_INDEX_COV2+ncovariate,false);
    belowthreshold=logreg.maximumLikelihoodRegression(param.iterations,param.threshold);

    cout<<"LR2 ->"<<markerid[markeridx]<<" ->";
    for (int i1=0; i1<logreg.beta.rows(); i1++)
      cout<<logreg.beta(i1)<<"; ";
    cout<<endl;

    results_value[RESULT_COLUMNS::STABLELRA]=belowthreshold?1:0;
    results_value[RESULT_COLUMNS::ORIO]=logreg.oddsratio(LR_INDEX_A1B0);
    results_value[RESULT_COLUMNS::ORIOL]=logreg.lowCI(LR_INDEX_A1B0);
    results_value[RESULT_COLUMNS::ORIOH]=logreg.highCI(LR_INDEX_A1B0);
    results_value[RESULT_COLUMNS::ORII]=logreg.oddsratio(LR_INDEX_A1B1);
    results_value[RESULT_COLUMNS::ORIIL]=logreg.lowCI(LR_INDEX_A1B1);
    results_value[RESULT_COLUMNS::ORIIH]=logreg.highCI(LR_INDEX_A1B1);
    results_value[RESULT_COLUMNS::OROI]=logreg.oddsratio(LR_INDEX_A0B1);
    results_value[RESULT_COLUMNS::OROIL]=logreg.lowCI(LR_INDEX_A0B1);
    results_value[RESULT_COLUMNS::OROIH]=logreg.highCI(LR_INDEX_A0B1);
    double apvalue=0;
    switch(param.apcalculation) {
      case PROPORTION_DISEASE:
        apvalue=logreg.calculateRERI(LR_INDEX_A1B0,LR_INDEX_A0B1,LR_INDEX_A1B1)/logreg.oddsratio(LR_INDEX_A1B1);
        break;
      case PROPORTION_EFFECT:
        apvalue=logreg.calculateRERI(LR_INDEX_A1B0,LR_INDEX_A0B1,LR_INDEX_A1B1)/logreg.oddsratio(LR_INDEX_A1B1)-1;
        break;
      case PROPORTION_CORRECTED:
        apvalue=logreg.calculateRERI(LR_INDEX_A1B0,LR_INDEX_A0B1,LR_INDEX_A1B1)/
                max(logreg.oddsratio(LR_INDEX_A1B1),logreg.oddsratio(LR_INDEX_A1B0)+logreg.oddsratio(LR_INDEX_A0B1)-1);
        break;
      }
    double aperror;
    aperror=logreg.APSEM(LR_INDEX_A1B0,LR_INDEX_A0B1,LR_INDEX_A1B1);
    results_value[RESULT_COLUMNS::AP]=apvalue;
    results_value[RESULT_COLUMNS::APL]=logreg.lowCI(apvalue,aperror);
    results_value[RESULT_COLUMNS::APH]=logreg.highCI(apvalue,aperror);
    boost::math::normal normaldist=boost::math::normal(0,aperror);
    results_value[RESULT_COLUMNS::APP]=1-cdf(normaldist,abs(apvalue))*2;
    }
  }
//------------------------------------------------------------------------------
void Analysis::alleleSummaryCount(int *alleles,int markeridx,int *phenotypex) {
  int offset;

  for (int i1=0; i1<nindividualid; i1++) {
    offset=(phenotypex[i1]==PHENOTYPE_AFFECTED?INDEX_CASE_PRIMARY:INDEX_CONTROL_PRIMARY);
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
void Analysis::setInteraction(int interactivemarkeridx) {
  double ratioriskalleleprimary,ratioriskallelesecondary;
  int riskhomozygote;
  int alleles[]={0,0,0,0};

  if (interactionfromfile!=NULL) {
    imarkinteraction=interactionfromfile;
    return;
    }
  if (imarkinteraction==NULL)
    imarkinteraction=new int[nindividualid];
  alleleSummaryCount(alleles,interactivemarkeridx,phenotype[ORIGINAL]);
  ratioriskalleleprimary=(double)alleles[INDEX_CASE_PRIMARY]/(double)alleles[INDEX_CONTROL_PRIMARY];
  ratioriskallelesecondary=(double)alleles[INDEX_CASE_SECONDARY]/(double)alleles[INDEX_CONTROL_SECONDARY];
  riskhomozygote=(ratioriskalleleprimary>ratioriskallelesecondary?HOMOZYGOTE_PRIMARY:HOMOZYGOTE_SECONDARY);
  for (int i1=0; i1<nindividualid; i1++) {
    imarkinteraction[i1]=NA_INTERACTION;
    if (!validGeneticData(i1,interactivemarkeridx,phenotype[ORIGINAL]))
      continue;
    switch (genotype[i1][interactivemarkeridx]) {
      case HOMOZYGOTE_PRIMARY:
        imarkinteraction[i1]=(riskhomozygote==HOMOZYGOTE_PRIMARY?INTERACTION:NO_INTERACTION);
        break;
      case HOMOZYGOTE_SECONDARY:
        imarkinteraction[i1]=(riskhomozygote==HOMOZYGOTE_SECONDARY?INTERACTION:NO_INTERACTION);
        break;
      case HETEROZYGOTE:
        imarkinteraction[i1]=(isDominantOrXMale(i1,interactivemarkeridx)?INTERACTION:NO_INTERACTION);
        break;
      }
    }
  }
//------------------------------------------------------------------------------
bool Analysis::validIndividualData(int individualidx,int markeridx,int *phenotypex) {
  return (genotype[individualidx][markeridx]!=ZYGOTE_UNKNOWN &&
          phenotypex[individualidx]!=PHENOTYPE_UNKNOWN &&
          interaction[individualidx]!=NA_INTERACTION);
  }
//------------------------------------------------------------------------------
bool Analysis::validGeneticData(int individualidx,int markeridx,int *phenotypex) {
  return (genotype[individualidx][markeridx]!=ZYGOTE_UNKNOWN &&
          phenotypex[individualidx]!=PHENOTYPE_UNKNOWN);
  }
//------------------------------------------------------------------------------
char Analysis::calculateRiskAllele(int markeridx, int *phenotypex, string *results) {
  char controlmaxallele,casemaxallele,caseminallele;
  int controlmax,casemax;
  int primarycount,secondarycount;
  double controlmaxratio,casemaxratio;
  int alleles[]={0,0,0,0};

  alleleSummaryCount(alleles,markeridx,phenotypex);
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
void Analysis::calculateRiskFactors(int markeridx,int * phenotypex,char riskallele,int recode) {
  int healthygenotype,riskgenotype;

  healthygenotype=riskallele==allele1[markeridx]?HOMOZYGOTE_SECONDARY:HOMOZYGOTE_PRIMARY;
  riskgenotype=riskallele==allele1[markeridx]?HOMOZYGOTE_PRIMARY:HOMOZYGOTE_SECONDARY;
  for (int i1=0; i1<nindividualid; i1++) {
    riskfactors[i1]=NA_RISK;
    if (!validIndividualData(i1,markeridx,phenotypex))
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
void Analysis::calculateRiskMatrix(int *phenotypex,int *riskmatrix) {
  int affstatus;

  fill_n(riskmatrix,N_RISK_MATRIX,0);
  for (int y1=0; y1<nindividualid; y1++) {
    fill_n(covariate1[y1],MATRIX_INDEX_COV1,NA_INTERACTION);
    fill_n(covariate2[y1],MATRIX_INDEX_COV2,NA_INTERACTION);
    if (riskfactors[y1]==NA_RISK)
      continue;
    fill_n(covariate1[y1],MATRIX_INDEX_COV1,NO_INTERACTION);
    fill_n(covariate2[y1],MATRIX_INDEX_COV2,NO_INTERACTION);
    affstatus=phenotypex[y1]-1;
    if (riskfactors[y1]==NO_RISK && interaction[y1]==NO_INTERACTION) {
      covariate2[y1][MATRIX_INDEX_A0B0]=INTERACTION;
      riskmatrix[affstatus]++;
      }
    if (riskfactors[y1]==RISK && interaction[y1]==NO_INTERACTION) {
      covariate2[y1][MATRIX_INDEX_A1B0]=INTERACTION;
      riskmatrix[RESULT_COLUMNS::IND10_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    if (riskfactors[y1]==NO_RISK && interaction[y1]>=INTERACTION) {
      covariate2[y1][MATRIX_INDEX_A0B1]=INTERACTION;
      riskmatrix[RESULT_COLUMNS::IND01_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    if (riskfactors[y1]==RISK && interaction[y1]>=INTERACTION) {
      covariate2[y1][MATRIX_INDEX_A1B1]=INTERACTION;
      riskmatrix[RESULT_COLUMNS::IND11_0-RESULT_COLUMNS::IND00_0+affstatus]++;
      }
    if (riskfactors[y1]==RISK)
      covariate1[y1][MATRIX_INDEX_A1m]=INTERACTION;
    if (interaction[y1]>=INTERACTION)
      covariate1[y1][MATRIX_INDEX_B1m]=INTERACTION;
    covariate1[y1][MATRIX_INDEX_A1mB1m]=covariate1[y1][MATRIX_INDEX_A1m]*covariate1[y1][MATRIX_INDEX_B1m];
    }
  }
//------------------------------------------------------------------------------
bool Analysis::belowCutOff(int *riskmatrix) {
  for (int i1=0; i1<N_RISK_MATRIX; i1++)
    if (riskmatrix[i1]<=param.cutoff)
      return true;
  return false;
  }
//------------------------------------------------------------------------------
void Analysis::setCleanData(int markeridx,int *y, double **x, VectorXd &desty, MatrixXd &destx, int dimx, bool a0b0) {
  int y2,x1,x2;

  desty.resize(nindividualid);
  destx.resize(nindividualid,dimx+1-!a0b0);
  for (int y1=y2=0; y1<nindividualid; y1++) {
    if (!validIndividualData(y1,markeridx,y))
      continue;
    desty(y2)=y[y1]-1;
    destx(y2,LR_BETA0)=1;
    for (x1=0,x2=1; x1<dimx; x1++) {
      if (x[y1][x1]==NA_INTERACTION)
        break;
      if (x1==MATRIX_INDEX_A0B0 && !a0b0)
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
  for (int i1=0; i1<nindividualid; i1++)
    if (interaction[i1]!=NA_INTERACTION)
      interaction[i1]=interaction[i1]==NO_INTERACTION?INTERACTION:NO_INTERACTION;
  }
//------------------------------------------------------------------------------


