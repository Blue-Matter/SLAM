
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

  // Annual time-series data (calculated from monthly)
  DATA_VECTOR(Effort_y_mean); // mean annual effort

  // Stock-recruit
  DATA_SCALAR(h); // steepness of BH-SRR

  // options
  DATA_INTEGER(Fit_Effort);
  DATA_INTEGER(Fit_CPUE);
  DATA_INTEGER(Fit_CAW);
  DATA_INTEGER(use_Eff_rwpen);  // use penalty for random walk in season effort
  DATA_INTEGER(use_R0rwpen);



  // ---- Estimated Parameters ----
  PARAMETER(ls50);  // log age-at-50% selectivity
  PARAMETER(lsdelta); // log interval age-50 - age-95% selectivity

  PARAMETER(logF_minit); // equilibrium fishing mortality for first age-classes

  PARAMETER(logq_effort);
  PARAMETER_VECTOR(logEffort_m_dev); // monthly mean fishing effort deviation (seasonal pattern)
  PARAMETER_VECTOR(logEffort_ts_dev); // fishing effort deviation for each timestep (month)
  PARAMETER(log_sigmaEff_m); // monthly rec dev sd (fixed or random effect)

  // random walk penalties for effort
  PARAMETER(log_Eff_m_SD); // standard deviation for random walk penalty for seasonal effort deviation

  PARAMETER_VECTOR(logR0_m_est); // average fraction of annual recruitment in each month
  PARAMETER(log_sigmaR0); // sd for random walk penalty for monthly recruitment

  PARAMETER_VECTOR(logRec_Devs); // monthly recruitment deviations
  PARAMETER(log_sigmaR); // monthly rec dev sd (fixed or random effect; usually fixed)

  // ---- Transform Parameters ----
  Type sigmaEff_m = exp(log_sigmaEff_m); // fishing effort monthly deviation sd
  Type sigmaR = exp(log_sigmaR); // rec process error dev sd
  Type sigmaR0 = exp(log_sigmaR0); // SD for random walk in R0_m (seasonal; monthly)

  Type Eff_m_SD = exp(log_Eff_m_SD); // standard deviation for random walk penalty for F

  // Selectivity-at-Age
  Type S50 = exp(ls50);
  Type Sdelta = exp(lsdelta);
  Type S95 = S50 + Sdelta;

  // ---- indexing variables ----
  int n_ages = Weight_Age.size(); // number of age classes
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

  // ---- Fishing effort and mortality ----
  // fishing mortality for initial age classes
  Type F_minit = exp(logF_minit);

  Type q_effort = exp(logq_effort);
  vector<Type> Effort_m_dev = exp(logEffort_m_dev);
  vector<Type> Effort_ts_dev = exp(logEffort_ts_dev);

  vector<Type> Effort_m(n_months); // fishing effort each timestep (month)
  vector<Type> F_m(n_months); // fishing mortality each timestep (month)
  Effort_m.setZero();
  F_m.setZero();
  int year_ind = -1;
  for(int m=0;m<n_months;m++){
    int m_ind = m % 12; // month index
    if (m_ind==0) year_ind = year_ind +1;
    Effort_m(m) = Effort_y_mean(year_ind) * Effort_m_dev(m_ind) * Effort_ts_dev(m);
    F_m(m) =  q_effort * Effort_m(m);
  }

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
  matrix<Type> B0_am(n_ages, 12); // B0 by age and month
  matrix<Type> SB0_am(n_ages, 12); // SB0 by age and month
  vector<Type> B0_m(12);
  vector<Type> SB0_m(12);
  B0_m.setZero();
  SB0_m.setZero();
  N_unfished.setZero();
  B0_am.setZero();
  SB0_am.setZero();

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
      B0_am(a, m_ind) =  N_unfished(a,m_ind) * Weight_Age(a) ;
      SB0_am(a, m_ind) =  N_unfished(a,m_ind) * Weight_Age(a) * Mat_at_Age(a);
    }
    B0_m(m_ind) = B0_am.col(m_ind).sum();
    SB0_m(m_ind) = SB0_am.col(m_ind).sum();
  }

  // ---- Equilibrium initial fished population ----
  matrix<Type> N_fished_eq(n_ages, 12);
  matrix<Type> SB_am_eq(n_ages, 12);
  vector<Type> SB_m_eq(12);
  N_fished_eq.setZero();
  SB_am_eq.setZero();
  SB_m_eq.setZero();

  vector<Type> Za_init(n_ages);
  vector<Type> Fa_init(n_ages);
  Fa_init.setZero(); // total mortality for initial age classes
  Za_init.setZero(); // total mortality for initial age classes
  for(int a=0;a<n_ages;a++){
    Fa_init(a) = F_minit * selA(a);
    Za_init(a) =  Fa_init(a) + M_at_Age(a);
  }

  for (int t=0; t<48; t++) {
    int m_ind = t % 12; // month index
    for(int a=1;a<n_ages;a++){

      if (t==0) {
        N_fished_eq(a,m_ind) = N_unfished(a-1,11) * exp(-Za_init(a-1)) * (1-PSM_at_Age(a-1));
      } else {
        if (m_ind==0) {
          N_fished_eq(a,m_ind) = N_fished_eq(a-1,11) * exp(-Za_init(a-1)) * (1-PSM_at_Age(a-1));
        } else {
          N_fished_eq(a,m_ind) = N_fished_eq(a-1,m_ind-1) * exp(-Za_init(a-1)) * (1-PSM_at_Age(a-1));
        }

      }
      SB_am_eq(a, m_ind) =  N_fished_eq(a,m_ind) * Weight_Age(a) * Mat_at_Age(a) * exp(-Fa_init(a)/2);;
    }
    SB_m_eq(m_ind) = SB_am_eq.col(m_ind).sum();
    N_fished_eq(0,m_ind) = BH_SRR(R0_m(m_ind), h, SB_m_eq(m_ind), SBpR);
  }

  // ---- Biomass matrices and vector ----
  matrix<Type> SB_am(n_ages, n_months); // SB by age and month
  SB_am.setZero();
  matrix<Type> B_am(n_ages, n_months); // B by age and month
  B_am.setZero();
  vector<Type> SB_m(n_months); // SB by month
  SB_m.setZero();
  vector<Type> B_m(n_months); // B by month
  B_m.setZero();
  matrix<Type> N_m(n_ages, n_months); // N by age and month
  N_m.setZero();

  for(int a=1;a<n_ages;a++){
    N_m(a,0) = N_fished_eq(a-1,11) * exp(-Za_init(a-1)) * (1-PSM_at_Age(a-1));
    SB_am(a,0) =  N_m(a,0) * Weight_Age(a) * Mat_at_Age(a)  * exp(-Fa_init(a)/2);
    B_am(a,0) = N_m(a,0) * Weight_Age(a) * exp(-F_ma(a,0)/2);
  }

  // recruitment in initial month
  SB_m(0) = SB_am.col(0).sum();
  B_m(0) = B_am.col(0).sum();
  N_m(0,0) = BH_SRR(R0_m(0), h, SB_m(0), SBpR) * exp(logRec_Devs(0) - pow(sigmaR,2)/Type(2.0));

  // ---- Population dynamics for remaining months ----
  for (int m=1; m<n_months; m++) {
    int m_ind = m % 12; // calendar month index
    for(int a=1;a<n_ages;a++){
      N_m(a,m) = N_m(a-1,m-1) * exp(-Z_ma(a-1, m-1)) * (1-PSM_at_Age(a-1));
      SB_am(a,m) = N_m(a,m) * Weight_Age(a) * Mat_at_Age(a) * exp(-F_ma(a,m)/2);
      B_am(a,m) = N_m(a,m) * Weight_Age(a) * exp(-F_ma(a,m)/2);
    }
    SB_m(m) = SB_am.col(m).sum();
    B_m(m) = B_am.col(m).sum();
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
  // mean 1 over time-steps where effort data exists
  Type Effmean = 0;
  Type Effsum = 0;
  Type Effn = 0;
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort(m)))) {
      Effsum += Effort_m(m);
      Effn += 1;
    }
    Effmean = Effsum/Effn;
  }

  // standardize to mean 1
  vector<Type> StEffort(n_months);
  StEffort.setZero();
  vector<Type> Effnll(n_months);
  Effnll.setZero();

  for (int m=0; m<n_months; m++) {
    StEffort(m) = Effort_m(m)/Effmean;
    if ((!R_IsNA(asDouble(Effort(m)))) & (Effort(m)!=0)) {
      Effnll(m)  -= dnorm_(log(StEffort(m)), log(Effort(m)), Effort_SD(m), true);
    }
  }

  // ---- Index ----
  vector<Type> predIndex(n_months);
  predIndex.setZero();

  // Calculate predicted Index
  predIndex = B_m;

  // mean 1 over time-steps where CPUE data exists
  Type CPUEmean = 0;
  Type CPUEsum = 0;
  Type CPUEn = 0;
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(CPUE(m)))) {
      CPUEsum += predIndex(m);
      CPUEn += 1;
    }
    CPUEmean = CPUEsum/CPUEn;
  }
  // standardize to mean 1
  vector<Type> CPUEnll(n_months);
  CPUEnll.setZero();
  vector<Type> stpredIndex(n_months);
  stpredIndex.setZero();
  for (int m=0; m<n_months; m++) {
    stpredIndex(m) = predIndex(m)/CPUEmean;
    if (!R_IsNA(asDouble(CPUE(m)))) {
      CPUEnll(m) -= dnorm_(log(stpredIndex(m)), log(CPUE(m)), CPUE_SD(m), true);
    }
  }

  // ---- Effort deviations ----
  Type effdevnll = 0;
  for(int m=0;m<n_months;m++){
    effdevnll -= dnorm(Effort_ts_dev(m), Type(0.0), sigmaEff_m, true);
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
  if (Fit_CAW>0) {
    nll_joint(0) =  CAWnll.sum();
  }
  // Effort
  if (Fit_Effort>0) {
    nll_joint(1) =  Effnll.sum();
  }

  // CPUE
  if (Fit_CPUE>0) {
    nll_joint(2) =  CPUEnll.sum();
  }

  // Effort deviations
  nll_joint(3) =  effdevnll;

  // Recruitment deviations
  nll_joint(4) =  recdevnll;

  // ---- Penalties ----

  // penalty for random walk in seasonal effort
  vector<Type> Eff_m_rwpen(12);
  Eff_m_rwpen.setZero();
  if (use_Eff_rwpen>0) {
    for (int m=1; m<12; m++) {
      Eff_m_rwpen(m-1) -= dnorm(Effort_m_dev(m), Effort_m_dev(m-1), Eff_m_SD, true);
    }
    Eff_m_rwpen(11) -= dnorm(Effort_m_dev(11), Effort_m_dev(0), Eff_m_SD, true);
  }
  nll_joint(5) =Eff_m_rwpen.sum();

  // penalty for random walk in logR0_m (seasonal recruitment)
  if (use_R0rwpen>0) {
    for(int m=1;m<12;m++){
      nll_joint(6) -= dnorm(logR0_m(m), logR0_m(m-1), sigmaR0, true);
    }
    nll_joint(6) -= dnorm(logR0_m(11), logR0_m(0), sigmaR0, true);
  }

  // ---- Total negative log-likelihood ----
  Type nll=0;
  nll = nll_joint.sum();

  // ---- Reports ----

  // Estimated parameters
  ADREPORT(F_m);
  ADREPORT(SPR);

  // Unfished
  REPORT(B0_am);
  REPORT(SB0_am);
  REPORT(B0_m);
  REPORT(SB0_m);
  REPORT(N_unfished);

  // Fished equilibirum
  REPORT(N_fished_eq);
  REPORT(SB_am_eq);
  REPORT(SB_m_eq);

  REPORT(F_minit);
  REPORT(Fa_init);
  REPORT(Za_init);

  REPORT(q_effort); // catchability

  // Predicted time-series
  REPORT(N_m); // numbers
  REPORT(SB_m); // spawning biomass
  REPORT(SB_am); // spawning biomass by age and month
  REPORT(predCB); // catch biomass
  REPORT(StEffort); // effort
  REPORT(stpredIndex); // Index of abundance
  REPORT(SPR); // SPR
  REPORT(Effort_m); // fishing effort
  REPORT(F_m); // fishing mortality
  REPORT(predCAW); // catch-at-weight

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
  REPORT(CPUEnll);
  REPORT(effdevnll);
  REPORT(recdevnll);
  REPORT(nll_joint);
  REPORT(nll);

  // Other stuff
  REPORT(AWK); // age-weight key
  REPORT(sigmaEff_m);
  REPORT(Effort_m_dev);
  REPORT(Effort_ts_dev);

  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
