#### build pomp model--test model ###

library(lubridate)
library(pomp)
library(doParallel)

statenames <- c("S","I","Vwp1","Vwp2","Vap","Vn","Awp","Anap","W","C")
t0 <- min(dates)

init.names<-c("S_0","I_0","Vwp1_0","Vwp2_0","Vap_0","Vn_0","Awp_0","Anap_0")
param.names<-c("Vwp1_wane_rate","Vwp2_wane_rate","Vap_wane_rate","Vn_wane_rate","rec_rate","Vwp_fail_rate","Vwp_symptom_rate","Vnap_fail_rate","Vnap_symptom_rate","beta0","beta1","beta_mod_Awp","beta_mod_Anap","rho","sigmaSE","lag","TDaP_waning_boost",init.names)

rproc <- Csnippet("
                  double N, beta, foi, dw, nbirths, Sbirths, Vwp_births1, Vwp_births2, Vap_births, wane_mod_TDaP;
                  double rate[24], trans[24];
                  
                  // pop size
                  N = S+I+Vwp1+Vwp2+Vap+Vn+Awp+Anap;
                  
                  // force of infection
                  beta = beta0 * (1 + beta1 * sin(2 * 3.14159 * (t-lag*365.25)/365));
                  foi = (beta*I+beta*beta_mod_Awp*Awp+beta*beta_mod_Anap*Anap)/N;
                  
                  // white noise (extra-demographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  // white noise (extra-demographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  // TDaP waning scalar
                  if (t<13149) {wane_mod_TDaP=1;}
                  if (t>=13149) {wane_mod_TDaP=TDaP_waning_boost;}
                  
                  // S
                  rate[0] = foi*dw/dt;  //infection rate (stochastic)
                  rate[1] = death_rate/365.25; // natural mortality
                  
                  // I
                  rate[2] = rec_rate;  // I recovery rate
                  rate[3] = death_rate/365.25; // natural mortality
                  
                  //Vwp1
                  rate[4] = Vwp_fail_rate*Vwp_symptom_rate*foi*dw/dt; // rate at which Vwp becomes I
                  rate[5] = Vwp_fail_rate*(1.0-Vwp_symptom_rate)*foi*dw/dt; // rate at which Vwp becomes Awp
                  rate[6] = Vwp1_wane_rate*wane_mod_TDaP;         // Vwp waning rate pre aP booster for 4th/5th dose
                  rate[7] = death_rate/365.25; // natural mortality
                  
                  //Vwp2
                  rate[8] = Vwp_fail_rate*Vwp_symptom_rate*foi*dw/dt; // rate at which Vwp becomes I
                  rate[9] = Vwp_fail_rate*(1.0-Vwp_symptom_rate)*foi*dw/dt; // rate at which Vwp becomes Awp
                  rate[10] = Vwp2_wane_rate*wane_mod_TDaP;         // Vwp waning rate pre aP booster for 4th/5th dose
                  rate[11] = death_rate/365.25; // natural mortality
                  
                  //Vap
                  rate[12] = Vnap_fail_rate*Vnap_symptom_rate*foi*dw/dt; // rate at which Vap becomes I
                  rate[13] = Vnap_fail_rate*(1.0-Vnap_symptom_rate)*foi*dw/dt; // rate at which Vap becomes Anap
                  rate[14] = Vap_wane_rate*wane_mod_TDaP;         // Vap waning rate pre Tdap booster for 4th/5th dose
                  rate[15] = death_rate/365.25; // natural mortality
                  
                  //Vn
                  rate[16] = Vnap_fail_rate*Vnap_symptom_rate*foi*dw/dt; // rate at which Vn becomes I
                  rate[17] = Vnap_fail_rate*(1.0-Vnap_symptom_rate)*foi*dw/dt; // rate at which Vn becomes Anap
                  rate[18] = Vn_wane_rate*wane_mod_TDaP;         // Vn waning rate
                  rate[19] = death_rate/365.25; // natural mortality
                  
                  //Awp
                  rate[20] = rec_rate;   // Awp recovery rate
                  rate[21] = death_rate/365.25; // natural mortality
                  
                  //Anap
                  rate[22] = rec_rate;   // Aap recovery rate
                  rate[23] = death_rate/365.25; // natural mortality
                  
                  // Poisson births
                  nbirths = rpois((birth_rate/365.25)*N*dt);
                  Sbirths = nearbyint((1.0-vacc_rate)*nbirths);
                  
                  // Vaccinations
                  if (t<7670) 
                  {
                  Vwp_births1 = nearbyint(vacc_rate*nbirths); 
                  Vwp_births2=0; 
                  Vap_births=0;
                  }
                  
                  if (t>=7670 && t<9862) 
                  {
                  Vwp_births1 = 0; 
                  Vwp_births2=nearbyint(vacc_rate*nbirths); 
                  Vap_births=0;
                  } 
                  
                  if (t>=9862) 
                  {
                  Vwp_births1 = 0; 
                  Vwp_births2=0; 
                  Vap_births=nearbyint(vacc_rate*nbirths);
                  } 
                  
                  // transitions between classes
                  reulermultinom(2, S, &rate[0], dt, &trans[0]);
                  reulermultinom(2, I, &rate[2], dt, &trans[2]);
                  reulermultinom(4, Vwp1, &rate[4], dt, &trans[4]);
                  reulermultinom(4, Vwp2, &rate[8], dt, &trans[8]);
                  reulermultinom(4, Vap, &rate[12], dt, &trans[12]);
                  reulermultinom(4, Vn, &rate[16], dt, &trans[16]);
                  reulermultinom(2, Awp, &rate[20], dt, &trans[20]);
                  reulermultinom(2, Anap, &rate[22], dt, &trans[22]);

                  S += Sbirths + trans[6] + trans[10] + trans[14] + trans[18] - trans[0] - trans[1];
                  I += trans[0] + trans[4] + trans[8] + trans[12] + trans[16] - trans[2] - trans[3];
                  Vwp1 += Vwp_births1 - trans[4] - trans[5] - trans[6] - trans[7];
                  Vwp2 += Vwp_births2 - trans[8] - trans[9] - trans[10] - trans[11];
                  Vap += Vap_births - trans[12] - trans[13] - trans[14] - trans[15];
                  Vn += trans[2] + trans[20] + trans[22] - trans[16] - trans[17] - trans[18] - trans[19];
                  Awp += trans[5] + trans[9] - trans[20] - trans[21];
                  Anap += trans[13] + trans[17] - trans[22] - trans[23];
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[0] + trans[4] + trans[8] + trans[12] + trans[16]; // true incidence
                  
                  if(I<0) {I=0;}
                  ")

rinit <- Csnippet("
                  double m = pop/(S_0+I_0+Vwp1_0+Vwp2_0+Vap_0+Vn_0+Awp_0+Anap_0);
                  S = nearbyint(m*S_0);
                  I = nearbyint(m*I_0);
                  Vwp1 = nearbyint(m*Vwp1_0);
                  Vwp2 = nearbyint(m*Vwp2_0);
                  Vap = nearbyint(m*Vap_0);
                  Vn = nearbyint(m*Vn_0);
                  Awp = nearbyint(m*Awp_0);
                  Anap = nearbyint(m*Anap_0);
                  W = 0.0;
                  C = 0.0;
                  ")

dmeas <- Csnippet("
                  lik = dpois(incidence,rho*C,give_log);
                  //if (!R_FINITE(lik)) {Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",S,I,Vwp1,Vwp2,Vap,Vn,Awp,Anap,beta0,beta1,lag,lik,C);}
                  ")

rmeas <- Csnippet("
                  incidence = rpois(rho*C);
                  if (incidence > 0.0) {
                  incidence = nearbyint(incidence);
                  } else {
                  incidence = 0.0;
                  }
                  ")


partrans=parameter_trans(
  log=c("TDaP_waning_boost","Vwp1_wane_rate","Vwp2_wane_rate","Vap_wane_rate","Vn_wane_rate","rec_rate","beta0","sigmaSE"),
  logit=c("Vwp_fail_rate","Vwp_symptom_rate","Vnap_fail_rate","Vnap_symptom_rate","beta_mod_Awp","beta_mod_Anap","rho","lag","beta1"),
  barycentric=c("S_0","I_0","Vwp1_0","Vwp2_0","Vap_0","Vn_0","Awp_0","Anap_0")
)


pomp(data=data,
     t0=min(dates),
     time="time",
     rprocess=euler(rproc,delta.t=1),
     dmeasure=dmeas,
     rmeasure=rmeas,
     rinit=rinit,
     partrans=partrans,
     covar=covariate_table(covartable,times="time"),
     accumvars=c("W","C"),
     statenames=statenames,
     paramnames=param.names,
     verbose = T
) -> m1

est.pars<-c("TDaP_waning_boost","Vwp1_wane_rate","Vwp2_wane_rate","Vap_wane_rate","Vn_wane_rate","rec_rate","beta0","sigmaSE",
            "Vwp_fail_rate","Vwp_symptom_rate","Vnap_fail_rate","Vnap_symptom_rate",
            "beta_mod_Awp","beta_mod_Anap","rho","lag","beta1","S_0","I_0","Vwp1_0","Vwp2_0","Vap_0","Vn_0","Awp_0","Anap_0")

rw.sd=rw.sd(TDaP_waning_boost=0.02,Vwp1_wane_rate=0.02,Vwp2_wane_rate=0.02,Vap_wane_rate=0.02,Vn_wane_rate=0.02,rec_rate=0.02,beta0=0.02,
            sigmaSE=0.02,Vwp_fail_rate=0.02,Vwp_symptom_rate=0.02,Vnap_fail_rate=0.02,Vnap_symptom_rate=0.02,
            beta_mod_Awp=0.02,beta_mod_Anap=0.02,
            rho=0.02,lag=0.02,beta1=0.02,S_0=0.02,I_0=0.02,Vwp1_0=0.02,Vwp2_0=0.02,Vap_0=0.02,Vn_0=0.02,Awp_0=0.02,Anap_0=0.02)

init.values.mat<-data.frame(S_0=LHS[,"S_0"],I_0=LHS[,"I_0"],Vwp1_0=LHS[,"Vwp1_0"],Vwp2_0=LHS[,"Vwp2_0"],Vap_0=LHS[,"Vap_0"],Vn_0=LHS[,"Vn_0"],Awp_0=LHS[,"Awp_0"],Anap_0=LHS[,"Anap_0"])

params.mat<-cbind(data.frame(TDaP_waning_boost = LHS[,"TDaP_waning_boost"],Vwp1_wane_rate = LHS[,"Vwp1_wane_rate"],Vwp2_wane_rate = LHS[,"Vwp2_wane_rate"],Vap_wane_rate = LHS[,"Vap_wane_rate"], Vn_wane_rate = LHS[,"Vn_wane_rate"], rec_rate = LHS[,"rec_rate"],
                             Vwp_fail_rate=LHS[,"Vwp_fail_rate"], Vwp_symptom_rate=LHS[,"Vwp_symptom_rate"], Vnap_fail_rate=LHS[,"Vnap_fail_rate"], Vnap_symptom_rate=LHS[,"Vnap_symptom_rate"], 
                             beta0 = LHS[,"beta0"], beta1 = LHS[,"beta1"],
                             beta_mod_Awp = LHS[,"beta_mod_Awp"], beta_mod_Anap = LHS[,"beta_mod_Anap"], rho=LHS[,"rho"],
                             sigmaSE=LHS[,"sigmaSE"],lag=LHS[,"lag"]),init.values.mat)

