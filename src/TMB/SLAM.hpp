
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

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
  DATA_VECTOR(EffExists); // logical (0 and 1) if effort data exists for this month
  DATA_INTEGER(nEffMonths); // total number of months

  DATA_VECTOR(CPUE); // monthly cpue - mean 1 over time-series
  DATA_VECTOR(CPUE_SD); // monthly cpue SD (log-space)

  // priors and penalties
  DATA_VECTOR(sigmaRprior); // prior on sigmaR
  DATA_VECTOR(F_meanprior); // prior on mean F

  // options
  DATA_INTEGER(Fit_Effort);
  DATA_INTEGER(Fit_CPUE);

  // Estimated Parameters (fixed)
  PARAMETER(log_sl50); // log length-at-50% selectivity
  PARAMETER(log_sldelta); // logSL95 - SL50

  PARAMETER_VECTOR(logR0_m_est); // monthly R0 - fraction
  PARAMETER(log_sigmaR0) // sd for random walk in monthly R0
  PARAMETER(logsigmaR); // monthly rec dev sd

  PARAMETER_VECTOR(logF_m); // monthly fishing mortality
  PARAMETER(log_sigmaF); // F standard deviation

  // Random Effects
  PARAMETER_VECTOR(logRec_Devs); // monthly recruitment deviations

  // Index variables
  int n_ages = Len_Age.size();
  int n_bins = LenMids.size();
  int n_months = CPUE.size();
  int ts_per_yr = 12.0;

  // Transform Parameters
  vector<Type> logR0_m(ts_per_yr);
  logR0_m.setZero();

  for(int m=1;m<12;m++){
    logR0_m(m) = logR0_m_est(m-1); // monthly mean rec
  }
  vector<Type> R0_m(ts_per_yr);
  R0_m.setZero();
  R0_m = exp(logR0_m);

  Type R0_mtotal = R0_m.sum();
  // standardize to sum to 1
  for(int m=0;m<ts_per_yr;m++){
    R0_m(m) = R0_m(m)/R0_mtotal;
  }

  vector<Type> F_m(n_months);
  F_m.setZero();
  F_m = exp(logF_m); // monthly fishing mortality

  Type sigmaF = exp(log_sigmaF); // F standard deviation
  Type sigmaR0 = exp(log_sigmaR0); // sd for monthly R0
  Type sigmaR = exp(logsigmaR); // monthly rec dev sd

  Type SL50 = exp(log_sl50);
  Type SLdelta = exp(log_sldelta);

  // Selectivity-at-Length
  vector<Type> selL(n_bins);
  selL.setZero();
  for(int l=0;l<n_bins;l++){
    selL(l) = 1 / (1 + exp(-log(Type(19))*(LenMids(l) - SL50)/SLdelta));
  }

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
  matrix<Type> N_m(n_ages, n_months);
  N_m.setZero();

  for(int a=0;a<n_ages;a++){
    if (a==0) {
      N_m(a,0) = R0_m(0)* exp(logRec_Devs(0) - pow(sigmaR,2)/Type(2.0));
    }
    if ((a>=1)) {
      N_m(a,0) = N_m(a-1,0) * exp(-Z_ma(a-1, 0)) * (1-PSM_at_Age(a-1));
    }
  }

  // loop over months
  for (int m=1; m<n_months; m++) {
    for(int a=0;a<n_ages;a++){
      if (a==0) {
        // month index
        int m_ind = m % 12;
        N_m(a,m) = R0_m(m_ind) * exp(logRec_Devs(m) - pow(sigmaR,2)/Type(2.0));
      }
      if ((a>=1)) {
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

  // Likelihood
  Type nll=0;
  vector<Type> nll_joint(8);
  nll_joint.setZero();


  // Effort
  vector<Type> RelEffort(n_months);
  RelEffort.setZero();
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort(m)))) {
      RelEffort(m) = F_m(m);
    }
  }

  // mean 1
  Type totEff = 0;
  totEff = RelEffort.sum();
  Type Effmean = 0;
  Effmean = totEff/nEffMonths;
  vector<Type> StEffort(n_months);
  StEffort.setZero();
  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort(m)))) {
      StEffort(m) = RelEffort(m)/Effmean;
    }
  }

  vector<Type> EffLike(n_months);
  EffLike.setZero();
  if (Fit_Effort>0) {
    for (int m=0; m<n_months; m++) {
      if (!R_IsNA(asDouble(Effort(m)))) {
        EffLike(m)  -= dnorm(log(StEffort(m)), log(Effort(m)), Effort_SD(m), true);
      }
    }
  }

  // CPUE
  vector<Type> predCPUE(n_months);
  predCPUE.setZero();
  if (Fit_CPUE>0) {
    for (int m=0; m<n_months; m++) {
      if (!R_IsNA(asDouble(Effort(m)))) {
        predCPUE(m) = predCB(m)/StEffort(m);
      }
    }
  }

  vector<Type> stpredCPUE(n_months);
  stpredCPUE.setZero();
  Type totCPUE = 0;
  totCPUE = predCPUE.sum();
  Type CPUEmean = 0;
  CPUEmean = totCPUE/nEffMonths;

  for (int m=0; m<n_months; m++) {
    if (!R_IsNA(asDouble(Effort(m)))) {
      stpredCPUE(m) = predCPUE(m)/CPUEmean;
    }
  }

  vector<Type> CPUELike(n_months);
  CPUELike.setZero();
  if (Fit_CPUE>0) {
    for (int m=0; m<n_months; m++) {
      if (!R_IsNA(asDouble(Effort(m)))) {
        CPUELike(m) -= dnorm(log(stpredCPUE(m)), log(CPUE(m)), CPUE_SD(m), true);
      }
    }
  }

  // CAL
  vector<Type> CALns(n_months);
  CALns.setZero();
  for (int m=0; m<n_months; m++) {
    CALns(m) = CAL.col(m).sum();
  }

  vector<Type> CALnll(n_months);
  CALnll.setZero();
  for (int m=0; m<n_months; m++) {
    if (CALns(m)>0) {
      vector<Type> prob(n_bins);
      vector<Type> CALm(n_bins);
      vector<Type> CALm_p(n_bins);
      vector<Type> CALm_pa(n_bins);

      prob.setZero();
      CALm.setZero();
      CALm_p.setZero();
      CALm_pa.setZero();

      prob = predCAL.col(m);
      CALm = CAL.col(m);
      CALm_p = CALm/CALm.sum();
      CALm_pa = CALm_p*CAL_ESS(m);

      CALnll(m) -= dmultinom(CALm_pa, prob, true);
    }
  }


  // likelihoods
  nll_joint(0) =  EffLike.sum();
  nll_joint(1) =  CALnll.sum();

  // rec devs
  for(int m=0;m<n_months;m++){
    nll_joint(2) -= dnorm(logRec_Devs(m), Type(0.0), sigmaR, true);
  }

  Type sigmaRpen;
  sigmaRpen = 0;
  sigmaRpen = Type(-1) * dnorm(log(sigmaR), log(sigmaRprior(0)), sigmaRprior(1), true);
  nll_joint(3) = sigmaRpen;


  // R0 pattern
  for(int m=1;m<ts_per_yr;m++){
    nll_joint(4) -= dnorm(logR0_m(m), logR0_m(m-1), sigmaR0, true);
  }
  nll_joint(4) -= dnorm(logR0_m(11), logR0_m(0), sigmaR0, true);

  // penalty random walk in F after initial
  for (int m=1; m<n_months; m++) {
    nll_joint(5) -= dnorm(F_m(m), F_m(m-1), sigmaF, true);
  }

  // penalty for mean F
  Type F_mean = 0;
  F_mean = F_m.sum()/F_m.size();
  nll_joint(6) = Type(-1)* dnorm(log(F_mean),log(F_meanprior(0)), F_meanprior(1), true);

  nll_joint(7) = CPUELike.sum();

  nll = nll_joint.sum();

  // Reports
  REPORT(SL50);
  REPORT(SLdelta);
  REPORT(F_m);
  REPORT(RelEffort);
  REPORT(predCAL);
  REPORT(sigmaF);
  REPORT(R0_m);
  REPORT(logRec_Devs);
  REPORT(sigmaR);
  REPORT(sigmaR0);
  REPORT(StEffort);
  REPORT(nll_joint);
  REPORT(N_m);
  REPORT(stpredCPUE);
  REPORT(CALnll);
  REPORT(Z_ma);
  REPORT(predCB);
  REPORT(selL);
  REPORT(selA);
  REPORT(ALK);
  REPORT(ALK_C);
  REPORT(predC_a);
  REPORT(EffLike);
  REPORT(CPUELike);

  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
