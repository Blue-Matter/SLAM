
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include "include/functions.hpp"

// Estimate optimal monthly fishing pattern given life-history,
// selectivity pattern, and monthly recruitment pattern

template<class Type>
Type optFpattern(objective_function<Type>* obj) {

  // Input information
  DATA_VECTOR(rec_pattern); // fraction of total annual recruitment in each month
  DATA_INTEGER(opt_type); // 0 to maximize yield, 1 to maximize HARA utility
  DATA_VECTOR(utilpow); // power for utility function

  // At-Age Schedules
  DATA_VECTOR(Len_Age);  // mean length at age
  DATA_VECTOR(SD_Len_Age); // standard deviation of length at age
  DATA_VECTOR(Wght_Age);  // mean weight at age
  DATA_VECTOR(Mat_at_Age);  // maturity at age
  DATA_VECTOR(M_at_Age); // natural mortality at age
  DATA_VECTOR(PSM_at_Age); // probability dying at-age (after spawning)

  DATA_VECTOR(LenBins);
  DATA_VECTOR(LenMids); // mid-points of the CAL bins

  // Estimated
  PARAMETER_VECTOR(logF_m); // log monthly fishing mortality

  int n_ages = Len_Age.size();
  int n_bins = LenMids.size();

  vector<Type> F_m(Type(12.0));
  F_m.setZero();
  F_m = exp(logF_m); // monthly fishing mortality

  // Selectivity-at-Length
  vector<Type> selL(n_bins);
  selL.setZero();
  for(int l=0;l<n_bins;l++){
    selL(l) = 1 / (1 + exp(-log(Type(19.0))*(LenMids(l) - SL50)/SLdelta));
  }

  // Generate Age-Length Key
  matrix<Type> ALK(n_ages, n_bins);
  ALK.setZero();
  ALK = generate_ALK(LenBins, Len_Age, SD_Len_Age, n_ages, n_bins);

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
  matrix<Type> F_ma(n_ages, Type(12.0));
  matrix<Type> M_ma(n_ages, Type(12.0));
  matrix<Type> Z_ma(n_ages, Type(12.0));
  F_ma.setZero();
  M_ma.setZero();
  Z_ma.setZero();

  for (int m=0; m<12; m++) {
    for(int a=0;a<n_ages;a++){
      F_ma(a,m) = F_m(m) * selA(a);
      M_ma(a,m) = M_at_Age(a);
      Z_ma(a,m) =  F_ma(a,m) +  M_ma(a,m);
    }
  }

  // numbers
  matrix<Type> N_m(n_ages, Type(36.0));  // run-out for 3 years to get rid of initial conditions
  N_m.setZero();

  for(int a=0;a<n_ages;a++){
    if (a==0) {
      N_m(a,0) = rec_pattern(0);
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
        N_m(a,m) = rec_pattern(m_ind);
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

  vector<Type> CB(Type(12.0)); // Catch biomass
  CB.setZero();
  CB = predCB.segment(25,36);

  Type nll=0;
  if (opt_type==0) {
    // maximize yield
    nll = -1*CB.sum();
  } else {
    // maximize HARA utility
    Type CBtotal = CB.sum();
    Type power = utilpow(0);
    Type HARA = pow(CBtotal, power);
    nll = -1*HARA
  }
  return(nll);

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
