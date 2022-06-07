
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type BHH_SRR(Type R0, Type h, Type SB, Type SBpR) {
  Type Rec = 0;
  Rec = (4*R0*h*SB)/(SBpR*R0*(1-h)+(5*h-1)*SB);
  return(Rec);
}


// Estimate optimal monthly fishing pattern given life-history,
// selectivity pattern, and monthly recruitment pattern

template<class Type>
Type optF(objective_function<Type>* obj) {

  // Input information
  DATA_VECTOR(rec_pattern); // fraction of total annual recruitment in each month
  DATA_INTEGER(opt_type); // 0 to maximize yield, 1 to maximize HARA utility
  DATA_VECTOR(utilpow); // power for utility function
  DATA_SCALAR(h); // steepness

  // At-Age Schedules
  DATA_VECTOR(Wght_Age);  // mean weight at age
  DATA_VECTOR(Mat_at_Age);  // maturity at age
  DATA_VECTOR(M_at_Age); // natural mortality at age
  DATA_VECTOR(PSM_at_Age); // probability dying at-age (after spawning)

  DATA_VECTOR(selA); // selectivity-at-age
  DATA_INTEGER(n_ages);


  // Estimated
  PARAMETER_VECTOR(logF_m); // log monthly fishing mortality

  vector<Type> F_m(12);
  F_m.setZero();
  F_m = exp(logF_m); // monthly fishing mortality

  // F, M, and Z by month and age
  matrix<Type> F_ma(n_ages, 12);
  matrix<Type> M_ma(n_ages, 12);
  matrix<Type> Z_ma(n_ages, 12);
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


  matrix<Type> N_m(n_ages, 12); // numbers
  N_m.setZero();
  vector<Type> SB(12);  // spawning biomass
  SB.setZero();
  vector<Type> Rec(12); // recruitment
  Rec.setZero();

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

  vector<Type> SB_a(n_ages);
  SB_a.setZero();
  // unfished - initialize
  for (int t=0; t<36; t++) { // run-out for 3 years to get rid of initial conditions
    int m_ind = t % 12; // month index
    for(int a=0;a<n_ages;a++){
      if (a==0) {
        N_m(a,m_ind) = rec_pattern(m_ind);
      } else {
        if (m_ind==0) {
          N_m(a,m_ind) = N_m(a-1,11) * exp(-M_ma(a-1, 11)) * (1-PSM_at_Age(a-1));
        } else {
          N_m(a,m_ind) = N_m(a-1,m_ind-1) * exp(-M_ma(a-1, m_ind-1)) * (1-PSM_at_Age(a-1));
        }
      }
      SB_a(a) = N_m(a,m_ind) * Wght_Age(a) * Mat_at_Age(a);
    }
    SB(m_ind) = SB_a.sum()*exp(-F_m(m_ind)/2);
    Rec(m_ind) = BHH_SRR(rec_pattern(m_ind), h, SB(m_ind), SBpR);
  }

  // with F and SRR
  SB_a.setZero();
  for (int t=0; t<60; t++) { // run out to ensure equilibrium
    int m_ind = t % 12; // month index
    for(int a=1;a<n_ages;a++){
      if (m_ind==0) {
        N_m(a,m_ind) = N_m(a-1,11) * exp(-Z_ma(a-1, 11)) * (1-PSM_at_Age(a-1));
      } else {
        N_m(a,m_ind) = N_m(a-1,m_ind-1) * exp(-Z_ma(a-1, m_ind-1)) * (1-PSM_at_Age(a-1));
      }
      SB_a(a) = N_m(a,m_ind) * Wght_Age(a) * Mat_at_Age(a);
    }
    SB(m_ind) = SB_a.sum()*exp(-F_m(m_ind)/2);
    Rec(m_ind) = BHH_SRR(rec_pattern(m_ind), h, SB(m_ind), SBpR);
    N_m(0,m_ind) = Rec(m_ind);
  }

  // Calculate catch
  matrix<Type> predC_a(n_ages, 12);
  matrix<Type> predCB_a(n_ages, 12);
  vector<Type> predCB(12);
  predC_a.setZero();
  predCB_a.setZero();
  predCB.setZero();
  for (int m=0; m<12; m++) {
    for(int a=0;a<n_ages;a++){
      predC_a(a,m) = N_m(a,m)*((1-Mat_at_Age(a))*exp(-M_ma(a,m)/2)+Mat_at_Age(a)*exp(-PSM_at_Age(a)/2))*(1-exp(-F_ma(a,m)));
      predCB_a(a,m) = predC_a(a,m) * Wght_Age(a);
    }
    predCB(m) = predCB_a.col(m).sum();
  }

  // vector<Type> CB(Type(12.0)); // Catch biomass
  // CB.setZero();
  // CB = predCB.segment(25,36);

  vector<Type> util(12);
  util.setZero();
  Type nll=0;
  Type CBtotal = predCB.sum();
  if (opt_type==0) {
    // maximize yield
    nll = -1*CBtotal;
  } else {
    // maximize HARA utility

    Type power = utilpow(0);
    for (int m=0; m<12; m++) {
      util(m) = pow(predCB(m), power);
    }
    Type HARA = util.sum();
    nll = -1*HARA;
  }

  REPORT(SBpR);
  REPORT(SB);
  REPORT(Rec);
  REPORT(util);
  REPORT(F_m);
  REPORT(predC_a);
  REPORT(predCB);
  REPORT(N_m);
  REPORT(Z_ma);
  REPORT(nll);
  return(nll);

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
