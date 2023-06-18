
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include "include/functions.h"

template<class Type>
Type SLAMv2(objective_function<Type>* obj) {

  // ---- Data ----
  // At-Age Schedules
  DATA_VECTOR(Weight_Age_Mean);  // mean weight at age
  DATA_VECTOR(Weight_Age_SD);  // standard deviation of weight at age (log-normal)
  DATA_VECTOR(Maturity_at_Age);  // maturity at age
  DATA_VECTOR(M_at_Age); // natural mortality at age
  DATA_VECTOR(Post_Spawning_Mortality); // probability dying at-age (after spawning)

  // Weight composition data
  DATA_MATRIX(CAW);    // CAW observations for each bin and month
  DATA_VECTOR(Weight_Bins);
  DATA_VECTOR(Weight_Mids); // mid-points of the CAW bins
  DATA_VECTOR(CAW_ESS); // number of independent observation of weight samples in each month

  // // Age composition data
  // DATA_MATRIX(CAA);    // CAA observations for each bin and month
  // DATA_VECTOR(CAA_ESS); // number of independent observation of age samples in each month

  // Monthly time-series data
  DATA_VECTOR(Effort_Mean); // monthly effort - mean 1 over time-series
  DATA_VECTOR(Effort_SD); // monthly effort SD (log-space)

  DATA_VECTOR(Index_Mean); // monthly cpue - mean 1 over time-series
  DATA_VECTOR(Index_SD); // monthly cpue SD (log-space)

  // Stock-recruit
  DATA_SCALAR(h); // steepness of BH-SRR

  // options
  DATA_INTEGER(Fit_Effort);
  DATA_INTEGER(Fit_Index);
  DATA_INTEGER(Fit_CAW);
  // DATA_INTEGER(Fit_CAA);
  DATA_INTEGER(use_Frwpen);
  DATA_INTEGER(use_R0rwpen);

  // ---- Estimated Parameters ----
  PARAMETER(ls50);  // log age-at-50% selectivity
  PARAMETER(lsdelta); // log interval age-50 - age-95% selectivity

  PARAMETER_VECTOR(logF_ts); // fishing mortality for each timestep (month)
  PARAMETER(log_sigmaF_m); // sd for random walk penalty for F

  PARAMETER(log_sigmaR0); // sd for random walk penalty for monthly recruitment

  PARAMETER_VECTOR(logRec_Devs); // monthly recruitment deviations
  PARAMETER(log_sigmaR); // monthly rec dev sd (fixed or random effect; usually fixed)

  // ---- Transform Parameters ----
  Type sigmaF_m = exp(log_sigmaF_m); // fishing effort monthly random walk sd
  Type sigmaR = exp(log_sigmaR); // rec process error dev sd
  Type sigmaR0 = exp(log_sigmaR0); // SD for random walk in R0_m (seasonal; monthly)

  // Selectivity-at-Age
  Type S50 = exp(ls50);
  Type Sdelta = exp(lsdelta);
  Type S95 = S50 + Sdelta;

  // ---- indexing variables ----
  int n_ages = Weight_Age_Mean.size(); // number of age classes
  int n_bins = Weight_Mids.size(); // number of size bins
  int n_months = Index_Mean.size(); // number of months of data

  // ---- Generate Age-Weight Key ----
  matrix<Type> AWK(n_ages, n_bins);
  AWK.setZero();
  AWK = generate_AWK(Weight_Bins, Weight_Age_Mean, Weight_Age_SD, n_ages, n_bins);


  // ---- Recruitment Deviations ----
  vector<Type> rec_devs(n_months); // recruitment deviations
  rec_devs.setZero();

  for (int m=0; m<n_months; m++) {
    rec_devs(m) = exp(logRec_Devs(m) - pow(sigmaR,2)/Type(2.0));
  }
  // ---- Seasonal Recruitment ----
  vector<Type> R0_m(12); // mean R0 for each calendar month
  R0_m.setZero();
  vector<Type> logR0_m(12); // mean R0 for each calendar month
  logR0_m.setZero();


  R0_m =  monthly_mean(rec_devs);
  Type R0_mtotal = R0_m.sum();

  // standardize to sum to 1
  for(int m=0;m<12;m++){
    R0_m(m) = R0_m(m)/R0_mtotal;
    logR0_m(m) = log(R0_m(m));
  }

  // ---- Fishing effort and mortality ----
  vector<Type> F_m(n_months); // fishing mortality each timestep (month)
  F_m.setZero();
  vector<Type> Effort_m(n_months); // relative predicted effort each timestep (month)
  Effort_m.setZero();

  for(int m=0;m<n_months;m++){
    F_m(m) =  exp(logF_ts(m));
    Effort_m(m) = F_m(m); // proportional to F
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
    surv0(a) = surv0(a-1)*exp(-M_ma(a-1,0))*(1-Post_Spawning_Mortality(a-1));
  }
  for (int a=0; a<n_ages; a++) {
    egg0(a) = surv0(a) * Weight_Age_Mean(a) * Maturity_at_Age(a);
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
        survF(a,m) = survF(a-1,m)*exp(-Z_ma(a-1, m)) * (1-Post_Spawning_Mortality(a-1));
      }
      eggFa(a) = survF(a,m) * Weight_Age_Mean(a) * Maturity_at_Age(a);
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
          N_unfished(a,m_ind) = N_unfished(a-1,11) * exp(-M_ma(a-1)) * (1-Post_Spawning_Mortality(a-1));
        } else {
          N_unfished(a,m_ind) = N_unfished(a-1,m_ind-1) * exp(-M_ma(a-1)) * (1-Post_Spawning_Mortality(a-1));
        }
      }
      B0_am(a, m_ind) =  N_unfished(a,m_ind) * Weight_Age_Mean(a) ;
      SB0_am(a, m_ind) =  N_unfished(a,m_ind) * Weight_Age_Mean(a) * Maturity_at_Age(a);
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

  // Calculate initial equilbrium F
  // equal to mean of estimated Fs
  Type F_minit = F_m.sum()/F_m.size();

  // add seasonal pattern for initial equilibrium Z
  vector<Type> Mean_monthly_F(12);
  Mean_monthly_F.setZero();

  Mean_monthly_F =  monthly_mean(F_m);

  // calculate average monthly pattern in F
  Type sumVal = Mean_monthly_F.sum();
  Type meanVal = sumVal/12;

  vector<Type> relMean_monthly_F(12);
  relMean_monthly_F.setZero();
  for (int m=0; m<12; m++) {
    relMean_monthly_F(m) =  Mean_monthly_F(m)/meanVal;
  }

  matrix<Type> Fa_init(n_ages, 12);
  Fa_init.setZero(); // total mortality for initial age classes
  matrix<Type> Za_init(n_ages, 12);
  Za_init.setZero(); // total mortality for initial age classes

  for (int m=0; m<12; m++) {
    for(int a=0;a<n_ages;a++){
      Fa_init(a, m) = F_minit * relMean_monthly_F(m) * selA(a);
      Za_init(a, m) =  Fa_init(a,m) + M_at_Age(a);
    }
  }


  for (int t=0; t<48; t++) {
    int m_ind = t % 12; // month index
    for(int a=1;a<n_ages;a++){
      if (t==0) {
        N_fished_eq(a,m_ind) = N_unfished(a-1,11) * exp(-Za_init(a-1, m_ind)) * (1-Post_Spawning_Mortality(a-1));
      } else {
        if (m_ind==0) {
          N_fished_eq(a,m_ind) = N_fished_eq(a-1,11) * exp(-Za_init(a-1, m_ind)) * (1-Post_Spawning_Mortality(a-1));
        } else {
          N_fished_eq(a,m_ind) = N_fished_eq(a-1,m_ind-1) * exp(-Za_init(a-1, m_ind)) * (1-Post_Spawning_Mortality(a-1));
        }
      }
      SB_am_eq(a, m_ind) =  N_fished_eq(a,m_ind) * Weight_Age_Mean(a) * Maturity_at_Age(a) * exp(-Fa_init(a, m_ind)/2);;
    }
    SB_m_eq(m_ind) = SB_am_eq.col(m_ind).sum();
    N_fished_eq(0,m_ind) = BH_SRR(R0_m(m_ind), h, SB_m_eq(m_ind), SBpR);
  }

  // ---- Biomass matrices and vector ----
  matrix<Type> SB_am(n_ages, n_months); // SB by age and month
  SB_am.setZero();
  matrix<Type> B_am(n_ages, n_months); // B by age and month
  B_am.setZero();
  matrix<Type> VB_am(n_ages, n_months); // VB by age and month
  VB_am.setZero();
  vector<Type> SB_m(n_months); // SB by month
  SB_m.setZero();
  vector<Type> B_m(n_months); // B by month
  B_m.setZero();
  vector<Type> VB_m(n_months); // VB by month
  VB_m.setZero();
  matrix<Type> N_m(n_ages, n_months); // N by age and month
  N_m.setZero();

  for(int a=1;a<n_ages;a++){
    N_m(a,0) = N_fished_eq(a-1,11) * exp(-Za_init(a-1)) * (1-Post_Spawning_Mortality(a-1));
    SB_am(a,0) =  N_m(a,0) * Weight_Age_Mean(a) * Maturity_at_Age(a)  * exp(-Fa_init(a)/2);
    B_am(a,0) = N_m(a,0) * Weight_Age_Mean(a);
    VB_am(a,0) = N_m(a,0) * Weight_Age_Mean(a) * selA(a);
  }

  // recruitment in initial month
  SB_m(0) = SB_am.col(0).sum();
  B_m(0) = B_am.col(0).sum();
  VB_m(0) = VB_am.col(0).sum();
  N_m(0,0) = BH_SRR(R0_m(0), h, SB_m(0), SBpR) * exp(logRec_Devs(0) - pow(sigmaR,2)/Type(2.0));

  // ---- Population dynamics for remaining months ----
  for (int m=1; m<n_months; m++) {
    int m_ind = m % 12; // calendar month index
    for(int a=1;a<n_ages;a++){
      N_m(a,m) = N_m(a-1,m-1) * exp(-Z_ma(a-1, m-1)) * (1-Post_Spawning_Mortality(a-1));
      SB_am(a,m) = N_m(a,m) * Weight_Age_Mean(a) * Maturity_at_Age(a) * exp(-F_ma(a,m)/2);
      B_am(a,m) = N_m(a,m) * Weight_Age_Mean(a);
      VB_am(a,m) = N_m(a,m) * Weight_Age_Mean(a) * selA(a);
    }
    SB_m(m) = SB_am.col(m).sum();
    B_m(m) = B_am.col(m).sum();
    VB_m(m) = VB_am.col(m).sum();
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
      predC_a(a,m) = N_m(a,m)*((1-Maturity_at_Age(a))*exp(-M_ma(a,m)/2)+Maturity_at_Age(a)*exp(-Post_Spawning_Mortality(a)/2))*(1-exp(-F_ma(a,m)));
      predCB_a(a,m) = predC_a(a,m) * Weight_Age_Mean(a);
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

  // Standardize predicted CAW to sum to one
  for (int m=0; m<n_months; m++) {
    Type temp = predCAW.col(m).sum();
    for(int w=0;w<n_bins;w++){
      predCAW(w,m) = predCAW(w,m)/temp;
    }
  }

  // Standardize predicted CAA to sum to one
  matrix<Type> predCAA(n_ages, n_months); // catch numbers by age and month
  predCAA.setZero();
  for (int m=0; m<n_months; m++) {
    Type temp2 = predC_a.col(m).sum();
    for(int a=0;a<n_ages;a++){
      predCAA(a,m) = predC_a(a,m)/temp2;
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

  // ---- Catch-at-Age ----
  // vector<Type> CAAns(n_months);
  // CAAns.setZero();
  // vector<Type> CAAnll(n_months);
  // CAAnll.setZero();
  //
  // for (int m=0; m<n_months; m++) {
  //   CAAns(m) = CAA.col(m).sum(); // sum of CAA observations
  //
  //   if (CAAns(m)>0) {
  //     // standardize observed CAA to sum 1
  //     vector<Type> CAAp_obs(n_ages);
  //     CAAp_obs.setZero();
  //     CAAp_obs = CAA.col(m)/CAA.col(m).sum();
  //
  //     // scale by effective sample size
  //     vector<Type> Ncaa_obs(n_ages);
  //     Ncaa_obs.setZero();
  //     Ncaa_obs = CAA_ESS(m) * CAAp_obs;
  //
  //     // multinomial likelihood
  //     vector<Type> predCAA_m(n_ages);
  //     predCAA_m.setZero();
  //     predCAA_m = predCAA.col(m);
  //     CAAnll(m) -= dmultinom_(Ncaa_obs, predCAA_m, true);
  //   }
  // }


  // ---- Relative Effort ----
  // mean 1 over time-steps where effort data exists
  Type Effmean = 0;
  Type Effsum = 0;
  Type Effn = 0;
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort_Mean(m)))) {
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
    if ((!R_IsNA(asDouble(Effort_Mean(m)))) & (Effort_Mean(m)!=0)) {
      Effnll(m)  -= dnorm_(log(StEffort(m)), log(Effort_Mean(m)), Effort_SD(m), true);
    }
  }

  // ---- Index ----
  vector<Type> predIndex(n_months);
  predIndex.setZero();

  // Calculate predicted Index
  predIndex = VB_m; // proportional to vulnerable biomass

  // mean 1 over time-steps where Index_Mean data exists
  Type CPUEmean = 0;
  Type CPUEsum = 0;
  Type CPUEn = 0;
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Index_Mean(m)))) {
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
    if (!R_IsNA(asDouble(Index_Mean(m)))) {
      CPUEnll(m) -= dnorm_(log(stpredIndex(m)), log(Index_Mean(m)), Index_SD(m), true);
    }
  }

  // ---- Recruitment deviations ----
  Type recdevnll = 0;
  for(int m=0;m<n_months;m++){
    recdevnll -= dnorm_(logRec_Devs(m), Type(0.0), sigmaR, true);
  }

  // ---- Joint likelihood ----
  vector<Type> nll_joint(7);
  nll_joint.setZero();

  // CAW
  if (Fit_CAW>0) {
    nll_joint(0) =  CAWnll.sum();
  }

  // // CAA
  // if (Fit_CAA>0) {
  //   nll_joint(1) =  CAAnll.sum();
  // }

  // Effort
  if (Fit_Effort>0) {
    nll_joint(1) =  Effnll.sum();
  }

  // Index
  if (Fit_Index>0) {
    nll_joint(2) =  CPUEnll.sum();
  }


  // Recruitment deviations
  if (sigmaR> 0.01) {
    nll_joint(3) =  recdevnll;
  }


  // ---- Penalties ----
  // penalty for random walk in F
  if (use_Frwpen>0) {
    for(int m=1;m<n_months;m++){
      nll_joint(4) -= dnorm_(logF_ts(m), logF_ts(m-1), sigmaF_m, true);
    }
  }

  // penalty for random walk in logR0_m (seasonal recruitment)
  if (use_R0rwpen>0) {
    for(int m=1;m<12;m++){
      nll_joint(5) -= dnorm_(logR0_m(m), logR0_m(m-1), sigmaR0, true);
    }
    nll_joint(5) -= dnorm_(logR0_m(11), logR0_m(0), sigmaR0, true);
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

  // Predicted time-series
  REPORT(N_m); // numbers
  REPORT(B_m); // biomass
  REPORT(SB_m); // spawning biomass
  REPORT(SB_am); // spawning biomass by age and month
  REPORT(predCB); // catch biomass
  REPORT(StEffort); // effort
  REPORT(stpredIndex); // Index of abundance
  REPORT(SPR); // SPR
  REPORT(F_m); // fishing mortality
  REPORT(predCAW); // catch-at-weight
  // REPORT(predCAA); // catch-at-age

  REPORT(predC_a);
  REPORT(predCB_a);
  REPORT(AWK);

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
  // REPORT(CAAnll);
  REPORT(Effnll);
  REPORT(CPUEnll);
  REPORT(recdevnll);
  REPORT(nll_joint);
  REPORT(nll);


  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
