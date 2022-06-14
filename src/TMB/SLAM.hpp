
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type calclog2(Type x) {
  Type ans = 0;
  ans = log(x)/log(2);
  return(ans);
}

template<class Type>
matrix<Type> generate_ALK(vector<Type> lbin, vector<Type> len_age, vector<Type> SD_LAA,
                          int n_age, int nlbin) {
  matrix<Type> ALK(n_age, nlbin);
  ALK.setZero();

  for(int a=0;a<n_age;a++) {
    for(int j=0;j<nlbin;j++) {
      if(j==nlbin-1) {
        ALK(a,j) = Type(1.0) - pnorm(lbin(j), len_age(a), SD_LAA(a));
      } else {
        ALK(a,j) = pnorm(lbin(j+1), len_age(a), SD_LAA(a));
        if(j>0) ALK(a,j) -= pnorm(lbin(j), len_age(a), SD_LAA(a));
      }
    }
  }
  return ALK;
}

template<class Type>
vector<Type> calSelL(vector<Type> Lens, Type LF5, Type LFS, Type Vmaxlen, Type Linf) {

  Type sigma_asc = (LFS-LF5)/sqrt(-calclog2(0.05));
  Type sigma_dec = 0;
  if (Vmaxlen<1) {
    sigma_dec= (Linf-LFS)/sqrt(-calclog2(Vmaxlen));
  } else {
    sigma_dec = 1E6;
  }


  int Lsize = Lens.size();
  vector<Type> sl(Lsize);
  sl.setZero();
  for (int i=0; i<Lsize; i++) {
    if (Lens(i)<=LFS) {
      Type temp = -pow((Lens(i)-LFS)/sigma_asc,2.0);
      sl(i) = pow(2, temp);
    } else{
      Type temp = -pow((Lens(i)-LFS)/sigma_dec,2.0);
      sl(i) = pow(2, temp);
    }
  }
  return(sl);

}


template<class Type>
Type SLAM(objective_function<Type>* obj) {

  // At-Age Schedules
  DATA_VECTOR(Len_Age);  // mean length at age
  DATA_VECTOR(SD_Len_Age); // standard deviation of length at age
  DATA_VECTOR(Wght_Age);  // mean weight at age
  DATA_VECTOR(Mat_at_Age);  // maturity at age
  DATA_VECTOR(M_at_Age); // natural mortality at age
  DATA_VECTOR(PSM_at_Age); // probability dying at-age (after spawning)

  // Catch-at-length data
  DATA_VECTOR(LenBins);
  DATA_VECTOR(LenMids); // mid-points of the CAL bins
  DATA_MATRIX(CAL);    // CAL observations for each bin and month
  DATA_VECTOR(CAL_ESS); // number of independent observation of length samples in each month

  // Monthly time-series data
  DATA_VECTOR(Effort); // monthly effort - mean 1 over time-series
  DATA_VECTOR(Effort_SD); // monthly effort SD (log-space)

  DATA_VECTOR(CPUE); // monthly cpue - mean 1 over time-series
  DATA_VECTOR(CPUE_SD); // monthly cpue SD (log-space)

  // priors and penalties
  DATA_VECTOR(F_meanprior); // prior on mean F

  // options
  DATA_INTEGER(Fit_Effort);
  DATA_INTEGER(Fit_CPUE);
  DATA_INTEGER(use_Fmeanprior);
  DATA_INTEGER(use_Frwpen);
  DATA_INTEGER(use_R0rwpen);

  DATA_SCALAR(Vmaxlen); // vulnerability at asymptotic length - fixed to 1 by default
  DATA_SCALAR(maxL); // length corresponding to Vmaxlen

  // Estimated Parameters
  PARAMETER(t_sl5); // log first length-at-5% selectivity
  PARAMETER(t_slfint); // log interval for first length-at-100% selectivity

  PARAMETER_VECTOR(logR0_m_est); // monthly R0 - fraction
  PARAMETER(log_sigmaR0); // sd for random walk in monthly R0

  PARAMETER_VECTOR(logF_m); // monthly fishing mortality
  PARAMETER(logF_minit); // mean fishing mortality for first age-classes

  PARAMETER_VECTOR(logRec_Devs); // monthly recruitment deviations

  //  variances
  PARAMETER(log_sigmaF); // F standard deviation
  PARAMETER(log_sigmaR); // monthly rec dev sd

  // Index variables
  int n_ages = Len_Age.size();
  int n_bins = LenMids.size();
  int n_months = CPUE.size();

  Type sigmaR = exp(log_sigmaR); // rec process error dev sd

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

  // Fishing mortality
  vector<Type> F_m(n_months);
  F_m.setZero();
  F_m = exp(logF_m); // monthly fishing mortality
  Type F_minit = exp(logF_minit); // mean fishing mortality in first time-step

  Type sigmaF = exp(log_sigmaF); // F standard deviation

  // Transform selectivity parameters
  Type SL5 = 0;
  Type SLFint = 0;
  SL5 = exp(t_sl5)/(1+exp(t_sl5)) * maxL;

  SLFint = exp(t_slfint)/(1+exp(t_slfint)) * maxL;
  Type SLFS = SL5 + SLFint;

  // Selectivity-at-Length
  vector<Type> selL(n_bins);
  selL.setZero();
  selL = calSelL(LenMids, SL5, SLFS, Vmaxlen, maxL);

  // Generate Age-Length Key
  matrix<Type> ALK(n_ages, n_bins);
  ALK.setZero();
  ALK = generate_ALK(LenBins, Len_Age, SD_Len_Age, n_ages, n_bins);

  // ALK for the catch
  matrix<Type> ALK_C(n_ages, n_bins);
  ALK_C.setZero();
  for(int a=0;a<n_ages;a++){
    for(int l=0;l<n_bins;l++){
      ALK_C(a,l) = ALK(a,l)*selL(l);
    }
  }

  for(int a=0;a<n_ages;a++){
    Type total = ALK_C.row(a).sum();
    ALK_C.row(a) = ALK_C.row(a)/total;
  }

  // Selectivity-at-Age
  vector<Type> selA(n_ages);
  selA.setZero();
  for(int a=0;a<n_ages;a++){
    vector<Type> temp(n_bins);
    temp.setZero();
    temp = ALK.row(a);
    for(int l=0;l<n_bins;l++){
      selA(a) += temp(l)*selL(l);
    }
  }

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

  // Population Dynamics - monthly time-step
  // first month
  matrix<Type> N_unfished(n_ages, 12);
  N_unfished.setZero();

  // initialize unfished population
  for (int t=0; t<36; t++) { // run-out for 3 years to get rid of initial conditions
    int m_ind = t % 12; // month index
    for(int a=0;a<n_ages;a++){
      if (a==0) {
        if (t==0) {
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

  for(int a=0;a<n_ages;a++){
    if (a==0) {
      N_m(a,0) = N_unfished(a, 0);
    }
    if (a>0) {
      N_m(a,0) = N_unfished(a-1,11) * exp(-Z_init(a-1)) * (1-PSM_at_Age(a-1));
    }
  }


  // loop over remaining months
  for (int m=1; m<n_months; m++) {
    int m_ind = m % 12; // calendar month index
    for(int a=0;a<n_ages;a++){
      if (a==0) {
        // month index
        N_m(a,m) = R0_m(m_ind) * exp(logRec_Devs(m) - pow(sigmaR,2)/Type(2.0));
      } else {
          N_m(a,m) = N_m(a-1,m-1) * exp(-Z_ma(a-1, m-1)) * (1-PSM_at_Age(a-1));
      }
    }
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

  // Calculate catch-at-length
  matrix<Type> predCAL(n_bins, n_months);
  predCAL.setZero();
  for (int m=0; m<n_months; m++) {
    for(int l=0;l<n_bins;l++){
      for(int a=0;a<n_ages;a++){
        predCAL(l,m) += predC_a(a,m)*ALK_C(a,l);
      }
    }
  }

  for (int m=0; m<n_months; m++) {
    Type temp = predCAL.col(m).sum();
    for(int l=0;l<n_bins;l++){
      predCAL(l,m) = predCAL(l,m)/temp;
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


  // Likelihoods

  // Catch-at-Length
  vector<Type> CALns(n_months);
  CALns.setZero();
  for (int m=0; m<n_months; m++) {
    CALns(m) = CAL.col(m).sum();
  }

  vector<Type> CALnll(n_months);
  CALnll.setZero();
  for (int m=0; m<n_months; m++) {
    if (CALns(m)>0) {
      vector<Type> CALp_obs(n_bins);
      CALp_obs.setZero();
      CALp_obs = CAL.col(m)/CAL.col(m).sum();
      N_obs = CAL_ESS(m) * CALp_obs;
      pred = predCAL.col(m);
      CALnll(m) -= (N_obs*log(pred)).sum();
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
      if ( stpredCPUE(m)>0) {
        CPUEnll(m) -= dnorm(log(stpredCPUE(m)), log(CPUE(m)), CPUE_SD(m), true);
      }
    }
  }

  // Recruitment deviations
  Type recdevnll = 0;
  for(int m=0;m<n_months;m++){
    recdevnll -= dnorm(logRec_Devs(m), Type(0.0), sigmaR, true);
  }

  vector<Type> nll_joint(8);
  nll_joint.setZero();

  nll_joint(0) =  CALnll.sum();
  if (Fit_Effort>0) {
    nll_joint(1) =  Effnll.sum();
  }
  if (Fit_CPUE>0) {
    nll_joint(2) =  CPUEnll.sum();
  }
  nll_joint(3) =  recdevnll;


  // Priors and penalties

  // penalty for mean F
  Type F_mean = 0;
  F_mean = F_m.sum()/F_m.size();
  if (use_Fmeanprior>0) {
    nll_joint(5) = Type(-1)* dnorm(log(F_mean),log(F_meanprior(0)), F_meanprior(1), true);
  }

  // penalty for random walk in F
  if (use_Frwpen>0) {
    for (int m=1; m<n_months; m++) {
      nll_joint(6) -= dnorm(F_m(m), F_m(m-1), sigmaF, true);
    }
  }

  // penalty for random walk in logR0_m
  if (use_R0rwpen>0) {
    for(int m=1;m<12;m++){
      nll_joint(7) -= dnorm(logR0_m(m), logR0_m(m-1), sigmaR0, true);
    }
    nll_joint(7) -= dnorm(logR0_m(11), logR0_m(0), sigmaR0, true);
  }

  Type nll=0;
  nll = nll_joint.sum();

  // Reports

  ADREPORT(SPR);


  REPORT(SPR);
  REPORT(ALK);
  REPORT(ALK_C);
  REPORT(SL5);
  REPORT(SLFS);
  REPORT(F_minit);
  REPORT(F_m);
  REPORT(R0_m);
  REPORT(sigmaR0);
  REPORT(logRec_Devs);
  REPORT(StEffort);
  REPORT(stpredCPUE);
  REPORT(predCAL);
  REPORT(predCB);
  REPORT(selL);
  REPORT(selA);

  REPORT(sigmaF);
  REPORT(sigmaR);

  REPORT(CALnll);
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
