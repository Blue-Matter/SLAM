
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
  DATA_VECTOR(Len_Age);
  DATA_VECTOR(SD_Len_Age);
  DATA_VECTOR(Wght_Age);
  DATA_VECTOR(Mat_at_Age);
  DATA_VECTOR(M_at_Age); // natural mortality
  DATA_VECTOR(PSM_at_Age); // probability dying at-age (after spawning)
  DATA_VECTOR(CAL_ESS);

  DATA_INTEGER(n_ages);
  // Monthly data
  DATA_VECTOR(CB); // catch-biomass by month

  DATA_INTEGER(n_months);
  DATA_SCALAR(sigma_C);

  // Assumed Known Parameters

  DATA_VECTOR(LenBins);
  DATA_VECTOR(LenMids);

  DATA_INTEGER(n_bins);
  DATA_INTEGER(ts_per_yr);

  DATA_MATRIX(CAL);

  // priors and penalties
  DATA_VECTOR(sigmaRprior);



  // Estimated Parameters (fixed)
  PARAMETER(log_sl50); // log length-at-50% selectivity
  PARAMETER(log_sldelta);


  PARAMETER_VECTOR(logR0_m); // monthly R0
  PARAMETER(logsigmaR); // monthly rec dev sd

  PARAMETER_VECTOR(logF_m); // monthly fishing mortality

  PARAMETER(log_sigmaF); // F standard deviation
  PARAMETER(log_sigmaR0) // sd for random walk in monthly R0

  // Random Effects
  PARAMETER_VECTOR(logRec_Devs); // monthly recruitment deviations

  // Transform Parameters
  vector<Type> R0_m(ts_per_yr);
  R0_m.setZero();
  R0_m = exp(logR0_m);

  vector<Type> F_m(n_months);
  F_m.setZero();
  F_m = exp(logF_m);

  Type sigmaF = exp(log_sigmaF);
  Type sigmaR0 = exp(log_sigmaR0);

  Type sigmaR = exp(logsigmaR);

  // Selectivity-at-Length
  Type SL50 = exp(log_sl50);
  Type SLdelta = exp(log_sldelta);
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
      Type total = ALK.row(a).sum();
      ALK_C(a,l) = ALK_C(a,l)/total;
    }
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

  // Population Dynamics - monthly time-step
  // first month
  matrix<Type> N_m(n_ages, n_months);
  N_m.setZero();

  for(int a=0;a<n_ages;a++){
    if (a==0) {
      N_m(a,0) = exp(R0_m(0)) * exp(logRec_Devs(0) - pow(sigmaR,2)/Type(2.0));
    }
    if ((a>=1)) {
      N_m(a,0) = N_m(a-1,0) * exp(-M_at_Age(a-1)+F_m(0)*selA(a-1)) * (1-PSM_at_Age(a-1));
    }
  }

  // loop over months
  for (int m=1; m<n_months; m++) {
    for(int a=0;a<n_ages;a++){
      if (a==0) {
        // month index
        int m_ind = m % 12;
        N_m(a,m) =  exp(R0_m(m_ind)) * exp(logRec_Devs(m) - pow(sigmaR,2)/Type(2.0));
      }
      if ((a>=1)) {
        N_m(a,m) = N_m(a-1,m-1) * exp(-M_at_Age(a-1)+F_m(m)*selA(a-1)) * (1-PSM_at_Age(a-1));
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
      predC_a(a,m) = N_m(a,m)*((1-Mat_at_Age(a))*exp(-M_at_Age(a)/2)+Mat_at_Age(a)*exp(-PSM_at_Age(a)/2))*(1-exp(-selA(a)*F_m(m)));
      predCB_a(a,m) = predC_a(a,m) * Wght_Age(a);
    }
    predCB(m) = predCB_a.col(m).sum();
  }

  vector<Type> logPC(n_months);
  logPC.setZero();
  // Catch likelihood
  for(int m=0;m<n_months;m++){
    if (CB(m)>0) {
      logPC(m) = dnorm(log(CB(m)), log(predCB(m)), sigma_C, true);
    }

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

  // likelihoods
  Type nll=0;
  vector<Type> nll_joint(6);
  nll_joint.setZero();

  nll_joint(0) = Type(-1) * logPC.sum();

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

  // CAL
  nll_joint(2) = CALnll.sum();

  // rec devs
  for(int m=0;m<n_months;m++){
    nll_joint(1) -= dnorm(logRec_Devs(m), Type(0.0), sigmaR, true);
  }

  Type sigmaRpen;
  sigmaRpen = 0;
  sigmaRpen = Type(-1) * dnorm(log(sigmaR), log(sigmaRprior(0)), sigmaRprior(1), true);
  nll_joint(3) = sigmaRpen;

  // penalty random walk in F after initial
  for (int m=2; m<n_months; m++) {
    nll_joint(4) -= dnorm(F_m(m), F_m(m-1), sigmaF, true);
  }

  // penalty random walk in monthly R0
  for (int m=1; m<ts_per_yr; m++) {
    nll_joint(5) -= dnorm(R0_m(m), R0_m(m-1), sigmaR0, true);
  }
  nll_joint(5) -= dnorm(R0_m(11), R0_m(0), sigmaR0, true);


  nll = nll_joint.sum();

  REPORT(SL50);
  REPORT(SLdelta);
  REPORT(selL);
  REPORT(predCB);
  REPORT(logRec_Devs);
  REPORT(N_m);
  REPORT(nll_joint);
  REPORT(R0_m);
  REPORT(F_m);
  REPORT(sigmaR);
  REPORT(predCAL);
  REPORT(CALnll);
  REPORT(sigmaF);
  REPORT(sigmaR0);


  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
