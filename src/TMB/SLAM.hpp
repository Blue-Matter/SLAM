
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include "include/functions.h"

template<class Type>
Type SLAM(objective_function<Type>* obj) {

  // ---- Data ----
  // At-Age Schedules
  DATA_VECTOR(Weight_Age);  // mean weight at age
  DATA_VECTOR(Weight_Age_SD);  // standard deviation of weight at age (log-normal)
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
  DATA_VECTOR(F_meanprior); // log-normal SD for monthly F

  // options
  DATA_INTEGER(Fit_Effort);
  DATA_INTEGER(Fit_CPUE);
  DATA_INTEGER(use_Frwpen);
  DATA_INTEGER(use_R0rwpen);
  DATA_INTEGER(use_Fmeanprior);

  // ---- Estimated Parameters ----
  PARAMETER(ls50);  // log age-at-50% selectivity
  PARAMETER(lsdelta); // log interval age-50 - age-95% selectivity

  PARAMETER_VECTOR(logF_minit); // equilibrium fishing mortality for first age-classes
  PARAMETER_VECTOR(logF_m); // monthly fishing mortality
  PARAMETER(log_sigmaF); // standard deviation for random walk penalty for F

  PARAMETER_VECTOR(logR0_m_est); // average fraction of annual recruitment in each month
  PARAMETER(log_sigmaR0); // sd for random walk penalty for monthly recruitment

  PARAMETER_VECTOR(logRec_Devs); // monthly recruitment deviations
  PARAMETER(log_sigmaR); // monthly rec dev sd (usually fixed)

  // ---- Transform Parameters ----
  Type sigmaR = exp(log_sigmaR); // rec process error dev sd
  Type sigmaR0 = exp(log_sigmaR0); // SD for random walk in R0_m
  Type sigmaF = exp(log_sigmaF); // standard deviation for random walk penalty for F

  // Selectivity-at-Age
  Type S50 = exp(ls50);
  Type Sdelta = exp(lsdelta);
  Type S95 = S50 + Sdelta;

  // ---- indexing variables ----
  int n_ages = Weight_Age.size(); // number of age classes
  int maxage = n_ages - 1;
  int n_bins = WghtMids.size(); // number of size bins
  int n_months = CPUE.size(); // number of months of data

  // ---- Generate Age-Weight Key ----
  matrix<Type> AWK(n_ages, n_bins);
  AWK.setZero();
  AWK = generate_AWK(WghtBins, Weight_Age, Weight_Age_SD, n_ages, n_bins);

  // ---- Seasonal Recruitment ----
  vector<Type> logR0_m(12); // R0 for each calendar month
  logR0_m.setZero();

  // map monthly mean rec fixing first month at log(1)
  for(int m=1;m<12;m++){
    logR0_m(m) = logR0_m_est(m-1);
  }

  vector<Type> R0_m(12);
  R0_m.setZero();
  R0_m = exp(logR0_m); // fraction of recruitment in each calendar month
  Type R0_mtotal = R0_m.sum();

  // standardize to sum to 1
  for(int m=0;m<12;m++){
    R0_m(m) = R0_m(m)/R0_mtotal;
  }

  // ---- Fishing mortality ----
  // fishing mortality for initial age classes
  vector<Type> F_minit(maxage);
  F_minit.setZero();
  F_minit = exp(logF_minit);

  // monthly fishing mortality
  vector<Type> F_m(n_months);
  F_m.setZero();
  F_m = exp(logF_m); // monthly fishing mortality

  // ---- Selectivity-at-Age ----
  vector<Type> selA(n_ages);
  selA.setZero();
  for(int a=0;a<n_ages;a++){
    selA(a) = 1 / (1 + exp(-log(Type(19))*((a - S50)/Sdelta)));
  }


  // --- Set up mortality matrices ----
  matrix<Type> M_ma(n_ages, n_months);
  matrix<Type> F_ma(n_ages, n_months);
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


  // ---- Calculate spawning biomass per recruit ----
  vector<Type> surv0(n_ages);
  surv0.setZero();
  surv0(0) = 1;
  vector<Type> egg0(n_ages);
  egg0.setZero();
  for (int a=1; a<n_ages; a++) {
    surv0(a) = surv0(a-1)*exp(-M_ma(a-1,0))*(1-PSM_at_Age(a-1));
  }
  for (int a=0; a<n_ages; a++) {
    egg0(a) = surv0(a) * Weight_Age(a) * Mat_at_Age(a);
  }
  Type SBpR = egg0.sum();

  // ---- Calculate SPR ----
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
      eggFa(a) = survF(a,m) * Weight_Age(a) * Mat_at_Age(a);
    }
    eggF(m) = eggFa.sum();
  }

  vector<Type> SPR(n_months);
  SPR.setZero();
  for (int m=0; m<n_months; m++) {
    SPR(m) = eggF(m)/SBpR;
  }


  // ---- Initialize unfished population ----
  // run-out for 3 years to get rid of initial conditions
  // and get equilibrium conditions for seasonal recruitment
  matrix<Type> N_unfished(n_ages, 12);
  N_unfished.setZero();
  for (int t=0; t<36; t++) {
    int m_ind = t % 12; // month index
    for(int a=0;a<n_ages;a++){
      if (a==0) {
        N_unfished(a,m_ind) = R0_m(m_ind);
      } else {
        if (m_ind==0) {
          N_unfished(a,m_ind) = N_unfished(a-1,11) * exp(-M_ma(a-1)) * (1-PSM_at_Age(a-1));
        } else {
          N_unfished(a,m_ind) = N_unfished(a-1,m_ind-1) * exp(-M_ma(a-1)) * (1-PSM_at_Age(a-1));
        }
      }
    }
  }

  // ---- Biomass matrices and vector ----
  matrix<Type> SB_am(n_ages, n_months); // SB by age and month
  SB_am.setZero();
  vector<Type> SB_m(n_months); // SB by month
  SB_m.setZero();

  // ---- Equilibrium initial fished population ----
  vector<Type> Z_init(maxage);
  Z_init.setZero(); // total mortality for initial age classes
  for(int a=0;a<maxage;a++){
    Z_init(a) =  F_minit(a) + M_at_Age(a);
  }

  matrix<Type> N_m(n_ages, n_months); // N by age and month
  N_m.setZero();

  for(int a=1;a<n_ages;a++){
    N_m(a,0) = N_unfished(a,11) * exp(-Z_init(a-1)) * (1-PSM_at_Age(a-1));
    SB_am(a,0) =  N_m(a,0) * Weight_Age(a) * Mat_at_Age(a); //  * exp(-F_minit(a)/2);
  }
  SB_m(0) = SB_am.col(0).sum(); // spawning biomass in initial month
  // recruitment in initial month
  N_m(0,0) = BH_SRR(R0_m(0), h, SB_m(0), SBpR) * exp(logRec_Devs(0) - pow(sigmaR,2)/Type(2.0));

  // ---- Population dynamics for remaining months ----
  for (int m=1; m<n_months; m++) {
    int m_ind = m % 12; // calendar month index
    for(int a=1;a<n_ages;a++){
      N_m(a,m) = N_m(a-1,m-1) * exp(-Z_ma(a-1, m-1)) * (1-PSM_at_Age(a-1));
      SB_am(a,m) = N_m(a,m) * Weight_Age(a) * Mat_at_Age(a) * exp(-F_ma(a,m)/2);
    }
    SB_m(m) = SB_am.col(m).sum();
    // recruitment
    N_m(0,m) = BH_SRR(R0_m(m_ind), h, SB_m(m), SBpR) * exp(logRec_Devs(m) - pow(sigmaR,2)/Type(2.0));
  }

  // ---- Calculate catch ----
  matrix<Type> predC_a(n_ages, n_months); // catch numbers by age and month
  matrix<Type> predCB_a(n_ages, n_months); // catch biomass by age and month
  vector<Type> predCB(n_months); // catch biomass by month
  predC_a.setZero();
  predCB_a.setZero();
  predCB.setZero();

  for (int m=0; m<n_months; m++) {
    for(int a=0;a<n_ages;a++){
      predC_a(a,m) = N_m(a,m)*((1-Mat_at_Age(a))*exp(-M_ma(a,m)/2)+Mat_at_Age(a)*exp(-PSM_at_Age(a)/2))*(1-exp(-F_ma(a,m)));
      predCB_a(a,m) = predC_a(a,m) * Weight_Age(a);
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

  // Standardize CAW data to sum to one
  for (int m=0; m<n_months; m++) {
    Type temp = predCAW.col(m).sum();
    for(int w=0;w<n_bins;w++){
      predCAW(w,m) = predCAW(w,m)/temp;
    }
  }

  // ---- Calculate likelihoods ----

  // --- Catch-at-Weight ---
  vector<Type> CAWns(n_months);
  CAWns.setZero();

  vector<Type> CAWnll(n_months);
  CAWnll.setZero();

  for (int m=0; m<n_months; m++) {
    CAWns(m) = CAW.col(m).sum(); // sum of CAW observations

    if (CAWns(m)>0) {
      // standardize observed CAW to sum 1
      vector<Type> CAWp_obs(n_bins);
      CAWp_obs.setZero();
      CAWp_obs = CAW.col(m)/CAW.col(m).sum();

      // scale by effective sample size
      vector<Type> N_obs(n_bins);
      N_obs.setZero();
      N_obs = CAW_ESS(m) * CAWp_obs;

      // multinomial likelihood
      vector<Type> pred(n_bins);
      pred.setZero();
      pred = predCAW.col(m);
      CAWnll(m) -= dmultinom_(N_obs, pred, true);
    }
  }

  // ---- Relative Effort ----
  // proportional to F
  vector<Type> RelEffort(n_months);
  RelEffort.setZero();
  for (int m=0; m<n_months; m++) {
      RelEffort(m) = F_m(m);
  }

  // mean 1 over time-steps where effort data exists
  Type Effmean = 0;
  Type Effsum = 0;
  Type Effn = 0;
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort(m)))) {
      Effsum += RelEffort(m);
      Effn += 1;
    }
    Effmean = Effsum/Effn;
  }

  // standardize to mean 1
  vector<Type> StEffort(n_months);
  StEffort.setZero();
  vector<Type> Effnll(n_months);
  Effnll.setZero();

  vector<Type> Effnll2(n_months);
  Effnll2.setZero();

  for (int m=0; m<n_months; m++) {
    StEffort(m) = RelEffort(m)/Effmean;
    if (!R_IsNA(asDouble(Effort(m)))) {
      Effnll(m)  -= dnorm(log(StEffort(m)), log(Effort(m)), Effort_SD(m), true);
      Effnll2(m)  -= dlognorm(StEffort(m), Effort(m), Effort_SD(m), true);
    }
  }

  // ---- CPUE ----
  vector<Type> predCPUE(n_months);
  predCPUE.setZero();

  // Calculate predicted CPUE
  for (int m=0; m<n_months; m++) {
      predCPUE(m) = predCB(m)/StEffort(m);
  }

  // mean 1 over time-steps where CPUE data exists
  Type CPUEmean = 0;
  Type CPUEsum = 0;
  Type CPUEn = 0;
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(CPUE(m)))) {
      CPUEsum += predCPUE(m);
      CPUEn += 1;
    }
    CPUEmean = CPUEsum/CPUEn;
  }
  // standardize to mean 1
  vector<Type> CPUEnll(n_months);
  CPUEnll.setZero();
  vector<Type> stpredCPUE(n_months);
  stpredCPUE.setZero();
  for (int m=0; m<n_months; m++) {
    stpredCPUE(m) = predCPUE(m)/CPUEmean;
    if (!R_IsNA(asDouble(CPUE(m)))) {
        CPUEnll(m) -= dnorm(log(stpredCPUE(m)), log(CPUE(m)), CPUE_SD(m), true);
    }
  }

  // ---- Recruitment deviations ----
  Type recdevnll = 0;
  for(int m=0;m<n_months;m++){
    recdevnll -= dnorm(logRec_Devs(m), Type(0.0), sigmaR, true);
  }

  // ---- Joint likelihood ----
  vector<Type> nll_joint(7);
  nll_joint.setZero();

  // CAW
  nll_joint(0) =  CAWnll.sum();

  // Effort
  if (Fit_Effort>0) {
    // nll_joint(1) =  Effnll.sum();
    nll_joint(1) =  Effnll2.sum();
  }

  // CPUE
  if (Fit_CPUE>0) {
    nll_joint(2) =  CPUEnll.sum();
  }

  // Recruitment deviations
  nll_joint(3) =  recdevnll;


  // ---- Priors and penalties ----

  // penalty for random walk in F
  vector<Type> Frwpen(n_months-1);
  Frwpen.setZero();
  if (use_Frwpen>0) {
    for (int m=1; m<n_months; m++) {
      Frwpen(m-1) -= dnorm(F_m(m), F_m(m-1), sigmaF, true);
    }
  }
  nll_joint(4) =Frwpen.sum();

  // penalty for random walk in logR0_m
  if (use_R0rwpen>0) {
    for(int m=1;m<12;m++){
      nll_joint(5) -= dnorm(logR0_m(m), logR0_m(m-1), sigmaR0, true);
    }
    nll_joint(5) -= dnorm(logR0_m(11), logR0_m(0), sigmaR0, true);
  }

  // penalty for mean monthly F
  vector<Type> F_month(12);
  F_month.setZero();
  vector<Type> F_count(12);
  F_count.setZero();
  // calculate mean F for months where CPUE data exists
  for (int m=0; m<n_months; m++) {
    int m_ind = m % 12; // calendar month index
    if (!R_IsNA(asDouble(CPUE(m)))) {
      F_month(m_ind) += F_m(m);
      F_count(m_ind) += Type(1);
    }
  }
  vector<Type> F_mean(12);
  F_mean.setZero();
  for (int m=0; m<12; m++) {
    F_mean(m) = F_month(m)/F_count(m);
  }

  vector<Type> F_month_NLL(12);
  F_month_NLL.setZero();
  if (use_Fmeanprior>0) {
    for (int m=0; m<n_months; m++) {
      int m_ind = m % 12; // calendar month index
      if (R_IsNA(asDouble(CPUE(m)))) {
        F_month_NLL(m_ind) += Type(-1)* dnorm(log(F_m(m)),log(F_mean(m_ind)), F_meanprior(0), true);
      }

    }
  }
  nll_joint(6) = F_month_NLL.sum();

  // ---- Total negative log-likelihood ----
  Type nll=0;
  nll = nll_joint.sum();


  // ---- Reports ----

  // Estimated parameters
  // ADREPORT(F_m);
  ADREPORT(SPR);

  // Predicted time-series
  REPORT(N_m); // numbers
  REPORT(SB_m); // spawning biomass
  REPORT(SB_am); // spawning biomass by age and month
  REPORT(predCB); // catch biomass
  REPORT(StEffort); // effort
  REPORT(stpredCPUE); // CPUE
  REPORT(SPR); // SPR
  REPORT(F_m); // fishing mortality
  REPORT(F_minit); // initial equilibrium Fs
  REPORT(predCAW); // catch-at-weight
  REPORT(N_unfished);

  // predicted seasonal recruitment
  REPORT(R0_m);
  REPORT(logRec_Devs); // recruitment deviations
  REPORT(sigmaR); // SD for rec devs


  // predicted selectivity-at-age
  REPORT(S50);
  REPORT(S95);
  REPORT(selA);

  // likelihoods
  REPORT(CAWnll);
  REPORT(Effnll);
  REPORT(Effnll2);
  REPORT(CPUEnll);
  REPORT(recdevnll);
  REPORT(nll_joint);
  REPORT(nll);

  // Other stuff
  REPORT(AWK); // age-weight key
  REPORT(sigmaF); // standard deviation for random walk penalty for F

  REPORT(Frwpen); // random walk F penalty
  REPORT(F_mean); // mean F
  REPORT(F_month_NLL); // mean F prior likelihood

  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
