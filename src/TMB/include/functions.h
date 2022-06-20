
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
matrix<Type> generate_AWK(vector<Type> WghtBins, vector<Type> Wght_Age,
                          vector<Type> Wght_Age_SD,
                          int n_age, int n_bins) {
  matrix<Type> AWK(n_age, n_bins);
  AWK.setZero();
  for(int a=0;a<n_age;a++) {
    Type mu = log(Wght_Age(a) -0.5*pow(Wght_Age_SD(a),2);
    for(int j=0;j<n_bins;j++) {
      if(j==n_bins-1) {
        AWK(a,j) = Type(1.0) - pnorm(log(WghtBins(j)),mu, Wght_Age_SD(a));
      } else {
        AWK(a,j) = pnorm(log(WghtBins(j+1)), mu, Wght_Age_SD(a));
        if(j>0) AWK(a,j) -= pnorm(log(WghtBins(j)), mu, Wght_Age_SD(a));
      }
    }
  }
  return AWK;
}

template<class Type>
Type BH_SRR(Type R0, Type h, Type SB, Type SBpR) {
  Type Rec = (4*R0*h*SB)/(SBpR*R0*(1-h)+(5*h-1)*SB);
  return Rec;
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
