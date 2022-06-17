
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include "include/functions.h"

template<class Type>
Type SLAM(objective_function<Type>* obj) {

  // At-Age Schedules
  DATA_VECTOR(Wght_Age);  // mean weight at age
  DATA_VECTOR(Wght_Age_SD);  // standard deviation of weight at age (log-normal)
  DATA_VECTOR(Mat_at_Age);  // maturity at age
  DATA_VECTOR(M_at_Age); // natural mortality at age
  DATA_VECTOR(PSM_at_Age); // probability dying at-age (after spawning)

  // Weight composition data
  DATA_VECTOR(WghtBins);
  DATA_VECTOR(WghtMids); // mid-points of the CAW bins
  DATA_MATRIX(CAW);    // CAL observations for each bin and month
  DATA_VECTOR(CAW_ESS); // number of independent observation of weight samples in each month

  // Monthly time-series data
  DATA_VECTOR(Effort); // monthly effort - mean 1 over time-series
  DATA_VECTOR(Effort_SD); // monthly effort SD (log-space)

  DATA_VECTOR(CPUE); // monthly cpue - mean 1 over time-series
  DATA_VECTOR(CPUE_SD); // monthly cpue SD (log-space)

  // Stock-recruit
  DATA_SCALAR(h); // steepness of BH-SRR

  // priors and penalties

  // options
  DATA_INTEGER(Fit_Effort);
  DATA_INTEGER(Fit_CPUE);
  DATA_INTEGER(use_Frwpen);
  DATA_INTEGER(use_R0rwpen);

  // Estimated Parameters
  PARAMETER(ls50);  //
  PARAMETER(lsdelta);

  PARAMETER_VECTOR(logR0_m_est); // monthly R0 - fraction
  PARAMETER(log_sigmaR0); // sd for random walk penalty in monthly R0

  PARAMETER_VECTOR(logF_m); // monthly fishing mortality
  PARAMETER(logF_minit); // mean fishing mortality for first age-classes

  PARAMETER_VECTOR(logRec_Devs); // monthly recruitment deviations

  //  variances
  PARAMETER(log_sigmaF); // standard deviation for random walk penalty for F
  PARAMETER(log_sigmaR); // monthly rec dev sd

  // Index variables
  int n_ages = Wght_Age.size();
  int n_bins = WghtMids.size();
  int n_months = CPUE.size();

  // Seasonal Recruitment
  vector<Type> logR0_m(12); // R0 for each calendar month
  logR0_m.setZero();
  for(int m=1;m<12;m++){
    logR0_m(m) = logR0_m_est(m-1); // map monthly mean rec
  }
  vector<Type> R0_m(12);
  R0_m.setZero();
  R0_m = exp(logR0_m);
  Type R0_mtotal = R0_m.sum();
  // standardize to sum to 1
  for(int m=0;m<12;m++){
    R0_m(m) = R0_m(m)/R0_mtotal;
  }

  Type sigmaR0 = exp(log_sigmaR0); // SD for random walk in R0_m
  Type sigmaR = exp(log_sigmaR); // rec process error dev sd

  // Fishing mortality
  vector<Type> F_m(n_months);
  F_m.setZero();
  F_m = exp(logF_m); // monthly fishing mortality
  Type F_minit = exp(logF_minit); // mean fishing mortality in first time-step

  Type sigmaF = exp(log_sigmaF); // F standard deviation

  // Transform selectivity parameters

  // Selectivity-at-Age
  Type S50 = exp(ls50);
  Type Sdelta = exp(lsdelta);
  Type S95 = S50 + Sdelta;

  vector<Type> selA(n_ages);
  selA.setZero();
  for(int a=0;a<n_ages;a++){
    selA(a) = 1 / (1 + exp(-log(Type(19))*(a - S50)/Sdelta));
  }

  // Generate Age-Weight Key
  matrix<Type> AWK(n_ages, n_bins);
  AWK.setZero();
  AWK = generate_AWK(WghtBins, Wght_Age, Wght_Age_SD, n_ages, n_bins);

  // F, M, and Z by month and age
  matrix<Type> F_ma(n_ages, n_months);
  matrix<Type> M_ma(n_ages, n_months);
  matrix<Type> Z_ma(n_ages, n_months);
  F_ma.setZero();
  M_ma.setZero();
  Z_ma.setZero();

  for (int m=0; m<n_months; m++) {
    for(int a=0;a<n_ages;a++){
      F_ma(a,m) = F_m(m) * selA(a);
      M_ma(a,m) = M_at_Age(a);
      Z_ma(a,m) =  F_ma(a,m) +  M_ma(a,m);
    }
  }

  // Calculate SPR
  // Calculate spawning biomass per recruit
  vector<Type> surv0(n_ages);
  surv0.setZero();
  surv0(0) = 1;
  vector<Type> egg0(n_ages);
  egg0.setZero();
  for (int a=1; a<n_ages; a++) {
    surv0(a) = surv0(a-1)*exp(-M_ma(a-1,0))*(1-PSM_at_Age(a-1));
  }
  for (int a=0; a<n_ages; a++) {
    egg0(a) = surv0(a) * Wght_Age(a) * Mat_at_Age(a);
  }
  Type SBpR = egg0.sum();


  matrix<Type> survF(n_ages, n_months);
  survF.setZero();
  vector<Type> eggFa(n_ages);
  eggFa.setZero();
  vector<Type> eggF(n_months);
  eggF.setZero();

  for (int m=0; m<n_months; m++) {
    for(int a=0;a<n_ages;a++){
      if (a==0) {
        survF(a,m) = 1;
      } else {
        survF(a,m) = survF(a-1,m)*exp(-Z_ma(a-1, m)) * (1-PSM_at_Age(a-1));
      }
      eggFa(a) = survF(a,m) * Wght_Age(a) * Mat_at_Age(a);
    }
    eggF(m) = eggFa.sum();
  }

  vector<Type> SPR(n_months);
  SPR.setZero();
  for (int m=0; m<n_months; m++) {
    SPR(m) = eggF(m)/SBpR;
  }

  // Population Dynamics - monthly time-step
  // first month
  matrix<Type> N_unfished(n_ages, 12);
  N_unfished.setZero();

  // initialize unfished population
  for (int t=0; t<36; t++) { // run-out for 3 years to get rid of initial conditions
    int m_ind = t % 12; // month index
    for(int a=0;a<n_ages;a++){
      if (a==0) {
        if (m_ind==0) {
          N_unfished(a,m_ind) = R0_m(m_ind);
        } else {
          N_unfished(a,m_ind) = R0_m(m_ind) * exp(logRec_Devs(m_ind) - pow(sigmaR,2)/Type(2.0));
        }
      } else {
        if (m_ind==0) {
          N_unfished(a,m_ind) = N_unfished(a-1,11) * exp(-M_ma(a-1, 11)) * (1-PSM_at_Age(a-1));
        } else {
          N_unfished(a,m_ind) = N_unfished(a-1,m_ind-1) * exp(-M_ma(a-1, m_ind-1)) * (1-PSM_at_Age(a-1));
        }
      }
    }
  }

  // first fished age-classes - month = 0
  // F, M, and Z by month and age
  vector<Type> F_init(n_ages);
  vector<Type> Z_init(n_ages);
  Z_init.setZero();
  F_init.setZero();

  for(int a=0;a<n_ages;a++){
    F_init(a) = F_minit * selA(a);
    Z_init(a) =  F_init(a) + M_at_Age(a);
  }
  matrix<Type> N_m(n_ages, n_months);
  N_m.setZero();

  matrix<Type> SB_am(n_ages, n_months);
  SB_am.setZero();

  vector<Type> SB_m(n_months);
  SB_m.setZero();

  for(int a=0;a<n_ages;a++){
    if (a==0) {
      N_m(a,0) = N_unfished(a, 0) * exp(logRec_Devs(0) - pow(sigmaR,2)/Type(2.0));
    }
    if (a>0) {
      N_m(a,0) = N_unfished(a-1,11) * exp(-Z_init(a-1)) * (1-PSM_at_Age(a-1));
      SB_am(a,0) = N_m(a,0) * Wght_Age(a) * Mat_at_Age(a) * exp(-F_init(a)/2);
    }
  }
  SB_m(0) = SB_am.col(0).sum();

  // loop over remaining months
  for (int m=1; m<n_months; m++) {
    int m_ind = m % 12; // calendar month index
    for(int a=1;a<n_ages;a++){
      N_m(a,m) = N_m(a-1,m-1) * exp(-Z_ma(a-1, m-1)) * (1-PSM_at_Age(a-1));
      SB_am(a,m) = N_m(a,m) * Wght_Age(a) * Mat_at_Age(a) * exp(-F_ma(a,m)/2);
    }
    SB_m(m) = SB_am.col(m).sum();
    // recruitment
    N_m(0,m) = BH_SRR(R0_m(m_ind), h, SB_m(m), SBpR) * exp(logRec_Devs(m) - pow(sigmaR,2)/Type(2.0));
  }

  // Calculate catch
  matrix<Type> predC_a(n_ages, n_months);
  matrix<Type> predCB_a(n_ages, n_months);
  vector<Type> predCB(n_months);
  predC_a.setZero();
  predCB_a.setZero();
  predCB.setZero();
  for (int m=0; m<n_months; m++) {
    for(int a=0;a<n_ages;a++){
      predC_a(a,m) = N_m(a,m)*((1-Mat_at_Age(a))*exp(-M_ma(a,m)/2)+Mat_at_Age(a)*exp(-PSM_at_Age(a)/2))*(1-exp(-F_ma(a,m)));
      predCB_a(a,m) = predC_a(a,m) * Wght_Age(a);
    }
    predCB(m) = predCB_a.col(m).sum();
  }

  // Calculate catch-at-weight
  matrix<Type> predCAW(n_bins, n_months);
  predCAW.setZero();
  for (int m=0; m<n_months; m++) {
    for(int w=0;w<n_bins;w++){
      for(int a=0;a<n_ages;a++){
        predCAW(w,m) += predC_a(a,m)*AWK(a,w);
      }
    }
  }

  for (int m=0; m<n_months; m++) {
    Type temp = predCAW.col(m).sum();
    for(int w=0;w<n_bins;w++){
      predCAW(w,m) = predCAW(w,m)/temp;
    }
  }



  // Likelihoods

  // Catch-at-Weight
  vector<Type> CAWns(n_months);
  CAWns.setZero();
  for (int m=0; m<n_months; m++) {
    CAWns(m) = CAW.col(m).sum();
  }

  vector<Type> CAWnll(n_months);
  CAWnll.setZero();
  for (int m=0; m<n_months; m++) {
    if (CAWns(m)>0) {
      vector<Type> CAWp_obs(n_bins);
      CAWp_obs.setZero();
      CAWp_obs = CAW.col(m)/CAW.col(m).sum();
      vector<Type> N_obs(n_bins);
      N_obs.setZero();
      N_obs = CAW_ESS(m) * CAWp_obs;
      vector<Type> pred(n_bins);
      pred.setZero();
      pred = predCAW.col(m);
      CAWnll(m) -= (N_obs*log(pred)).sum();
    }
  }

  // Effort
  vector<Type> RelEffort(n_months);
  RelEffort.setZero();
  Type nEffMonths = 0;
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort(m)))) {
      RelEffort(m) = F_m(m);
      nEffMonths += 1;
    }
  }

  // mean 1
  Type Effmean = 0;
  Effmean = RelEffort.sum()/nEffMonths;
  vector<Type> StEffort(n_months);
  StEffort.setZero();
  vector<Type> Effnll(n_months);
  Effnll.setZero();
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort(m)))) {
      StEffort(m) = RelEffort(m)/Effmean;
      Effnll(m)  -= dnorm(log(StEffort(m)), log(Effort(m)), Effort_SD(m), true);
    }
  }

  // CPUE
  vector<Type> predCPUE(n_months);
  predCPUE.setZero();
  vector<Type> stpredCPUE(n_months);
  stpredCPUE.setZero();
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort(m)))) {
      predCPUE(m) = predCB(m)/StEffort(m);
    }
  }

  Type CPUEmean = 0;
  CPUEmean = predCPUE.sum()/nEffMonths;
  vector<Type> CPUEnll(n_months);
  CPUEnll.setZero();
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort(m)))) {
      stpredCPUE(m) = predCPUE(m)/CPUEmean;
      if (stpredCPUE(m)>0) {
        CPUEnll(m) -= dnorm(log(stpredCPUE(m)), log(CPUE(m)), CPUE_SD(m), true);
      }
    }
  }

  // Recruitment deviations
  Type recdevnll = 0;
  for(int m=0;m<n_months;m++){
    recdevnll -= dnorm(logRec_Devs(m), Type(0.0), sigmaR, true);
  }

  vector<Type> nll_joint(6);
  nll_joint.setZero();

  nll_joint(0) =  CAWnll.sum();
  if (Fit_Effort>0) {
    nll_joint(1) =  Effnll.sum();
  }
  if (Fit_CPUE>0) {
    nll_joint(2) =  CPUEnll.sum();
  }
  nll_joint(3) =  recdevnll;


  // Priors and penalties

  // penalty for random walk in F
  if (use_Frwpen>0) {
    for (int m=1; m<n_months; m++) {
      nll_joint(4) -= dnorm(F_m(m), F_m(m-1), sigmaF, true);
    }
  }

  // penalty for random walk in logR0_m
  if (use_R0rwpen>0) {
    for(int m=1;m<12;m++){
      nll_joint(5) -= dnorm(logR0_m(m), logR0_m(m-1), sigmaR0, true);
    }
    nll_joint(5) -= dnorm(logR0_m(11), logR0_m(0), sigmaR0, true);
  }

  Type nll=0;
  nll = nll_joint.sum();

  // Reports

  ADREPORT(SPR);
  ADREPORT(S50);
  ADREPORT(S95);

  REPORT(SPR);
  REPORT(AWK);

  REPORT(S50);
  REPORT(S95);

  REPORT(SB_m);
  REPORT(F_minit);
  REPORT(F_m);
  REPORT(R0_m);
  REPORT(sigmaR0);
  REPORT(logRec_Devs);
  REPORT(StEffort);
  REPORT(stpredCPUE);
  REPORT(predCAW);
  REPORT(predCB);
  REPORT(selA);

  REPORT(sigmaF);
  REPORT(sigmaR);

  REPORT(CAWnll);
  REPORT(Effnll);
  REPORT(CPUEnll);
  REPORT(recdevnll);
  REPORT(nll_joint);

  REPORT(N_unfished);
  REPORT(N_m);

  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
