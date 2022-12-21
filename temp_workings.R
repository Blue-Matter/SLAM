
devtools::install_github('Blue-Matter/SLAM',
                         auth_token='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')


template<class Type>
  matrix<Type> generate_AWK(vector<Type> WghtBins, vector<Type> Wght_Age,
                            vector<Type> Wght_Age_SD,
                            int n_age, int n_bins) {
    matrix<Type> AWK(n_age, n_bins);
    AWK.setZero();

    for(int a=0;a<n_age;a++) {
      for(int j=0;j<n_bins;j++) {
        if(j==n_bins-1) {
          ALK(a,j) = Type(1.0) - pnorm(lbin(j), len_age(a), SD_LAA(a));
        } else {
          ALK(a,j) = pnorm(lbin(j+1), len_age(a), SD_LAA(a));
          if(j>0) ALK(a,j) -= pnorm(lbin(j), len_age(a), SD_LAA(a));
        }
      }
    }
    return ALK;
  }


library(SLAM)
Pars <- Load_Scenario(1)
WghtCV <- 0.24
Ages <- 0:Pars$maxage
nAge <- length(Ages)
LatAge <- Pars$Linf*(1-exp(-Pars$K*(Ages-Pars$t0)))
LatAge[LatAge==0] <- 0.001
LenSD <- LatAge * Pars$LenCV
Pars$wa <- 1E-5
WatAge <- Pars$wa*LatAge^Pars$wb
WatAge_SD <- rep(WghtCV, length(WatAge))

Pars$BinWidth <- 0.5
Pars$BinMax <- 15
# Probability of weight-at-age
Wght_Bins <- seq(0, to=Pars$BinMax, by=Pars$BinWidth)
nbins <- length(Wght_Bins)-1
By <- Wght_Bins[2]-Wght_Bins[1]
Wght_Mids <- seq(Wght_Bins[1]+0.5*By, by=By, length.out=nbins)

AWK <- matrix(0, nrow=nAge, ncol=nbins)
mu <- log(WatAge) -0.5*WatAge_SD^2

pnorm((log(Wght_Bins[2])-log(WatAge))/WatAge_SD, 0, 1)
pnorm((log(Wght_Bins[2])-mu)/WatAge_SD, 0, 1) == plnorm(Wght_Bins[2], mu, WatAge_SD)


AWK[,1] <- plnorm(Wght_Bins[2], mu, WatAge_SD) # probability of weight-at-age
for (i in 2:(nbins-1)) {
  AWK[,i] <- plnorm(Wght_Bins[i+1], mu, WatAge_SD) -
    plnorm(Wght_Bins[i], mu, WatAge_SD)
}
AWK[,nbins] <- 1 - plnorm(Wght_Bins[nbins], mu, WatAge_SD)


AWK2 <- matrix(0, nrow=nAge, ncol=nbins)
AWK2[,1]
mu <- log(WatAge)  -0.5*WatAge_SD^2
pnorm((log(Wght_Bins[2])-log(WatAge))/WatAge_SD, 0, 1) # probability of weight-at-age
for (i in 2:(nbins-1)) {
  AWK2[,i] <- pnorm(log(Wght_Bins[i+1]), mu, WatAge_SD) -
    pnorm(log(Wght_Bins[i]), mu, WatAge_SD)
}
AWK2[,nbins] <- 1 - pnorm(log(Wght_Bins[nbins]), mu, WatAge_SD)


a <- 10
plot(Wght_Mids, AWK[a,], type='b', pch=16)
lines(Wght_Mids, AWK2[a,], col='red')

tt <- exp(rnorm(1E6, log(WatAge[a])-0.5*WatAge_SD^2 , WatAge_SD))
Wght_Bins2 <- Wght_Bins
Wght_Bins2[length(Wght_Bins2)] <- max(tt)*1.1
t1 <- table(cut(tt, Wght_Bins2))
t1 <- t1/sum(t1)
lines(Wght_Mids, t1, col='blue', type='b', pch=16)

mean(tt)
WatAge[a]
sd(log(tt))



t1<- plnorm(Wght_Mids[15], log(WatAge), WatAge_SD) - plnorm(Wght_Mids[14], log(WatAge), WatAge_SD)
t2 <- pnorm(log(Wght_Mids[15]), log(WatAge), WatAge_SD) - pnorm(log(Wght_Mids[14]), log(WatAge), WatAge_SD)


a <- 3
tt <- dlnorm(Wght_Mids, log(WatAge[a]), WatAge_SD[a])
plot(Wght_Mids, tt, type='l')

t1 <- plnorm(Wght_Mids[2], log(WatAge), WatAge_SD)
plot(Wght_Mids, t1, type='l')

plot(WatAge, type='l')


tt <- plnorm((log(Wght_Bins[a]) - log(WatAge))/WatAge_SD, 0, 1)
plot(tt, type='l')












# ---- Compare Expected Yield from Assessment vs Operating Model ----

scen_list <- list()
for (scen in 1:nrow(grid)) {
  message('Scenario: ', scen)
  sub_grid <- grid[scen,]
  Name <- paste0(apply(sub_grid, 2, as.character), collapse="_")
  fl <-  file.path(paste0('Results/Sim_Test/', Name, '.rda'))
  obj <- readRDS(fl)

  sim_list<-list()
  for (sim in 1:nsim) {
    message("Sim: ", sim)
    SimPop <- obj[[sim]]$Sim
    Assess_obj <- obj[[sim]]$Assess

    # True
    Data_true <- list()
    Data_true$Weight_Age <- SimPop$Pars$Weight_Age
    Data_true$R0 <- SimPop$Pars$Rbar
    nage <- length(Data_true$Weight_Age)
    Data_true$Mat_at_Age <- SimPop$Pars$Mat_at_Age
    Data_true$M_at_Age <- rep(SimPop$Pars$M, nage)
    Data_true$PSM_at_Age <- SimPop$Pars$PSM_at_Age

    # Run simulation model with this pattern in F and compare yield
    opt1 <- Optimize(Data_true, SimPop$Rec_Pattern, SimPop$Sel_at_Age, opt_type=1,
                     utilpow=1,
                     assumed_h=SimPop$Pars$h)

    # Assessment
    opt2 <- Optimize(Data_true, Assess_obj$rep$R0_m,
                     Assess_obj$rep$selA,
                     opt_type=1,
                     utilpow=1,
                     assumed_h=SimPop$Pars$h)

    # Long-term expected yield fishing at assessment's recommendation
    Pars <- SimPop$Pars
    Pars$rec_devs[,1] <- 1
    Pars$Effort$Effort <- opt1$F_m
    Pars$q <- 1
    sim1 <- Simulate(Pars)

    Pars <- SimPop$Pars
    Pars$rec_devs[,1] <- 1
    Pars$Effort$Effort <- opt2$F_m
    Pars$q <- 1
    sim2 <- Simulate(Pars)

    sum(opt1$predCB*Pars$Rbar)
    sum(opt2$predCB*Pars$Rbar)
    sum(sim1$Catch_Biomass[(nts-11):nts])
    sum(sim2$Catch_Biomass[(nts-11):nts])


    sim_list[[sim]] <- data.frame(OM=sum(opt1$predCB),
                                  Est=sum(opt2$predCB), Sim=sim,
                                  Rec_Scen=sub_grid$rec_scen_names,
                                  RecSD=sub_grid$recsd_scen_names,
                                  M=sub_grid$M_scen_names)
  }
  scen_list[[scen]] <- do.call('rbind', sim_list)

}
DF_compare <- do.call('rbind', scen_list)

DF_compare$RE <- (DF_compare$Est-DF_compare$OM)/DF_compare$OM


ggplot(DF_compare, aes(x=Rec_Scen, y=RE, fill=as.factor(M))) +
  facet_grid(~RecSD) +
  geom_boxplot()

DF_compare %>% filter(Rec_Scen=='Constant', M==0.1, RecSD=="High")

DF_compare %>% filter(Rec_Scen=='Pulse', M==0.1, RecSD=="High", Sim==1)
DF_compare %>% filter(Rec_Scen=='Diffuse', M==0.1, RecSD=="High", Sim==1)

