#### build pomp model--test model ###

library(lubridate)
library(pomp)
library(doParallel)

statenames <- c("S","I","Vwp","Vap","Vn","Awp","Aap","An","W","C")
t0 <- min(dates)

init.names<-c("S_0","I_0","Vwp_0","Vap_0","Vn_0","Awp_0","Aap_0","An_0")
param.names<-c("Vwp_wane_rate","Vap_wane_rate","Vn_wane_rate","rec_rate","Vwp_fail_rate","Vwp_symptom_rate","Vap_fail_rate","Vap_symptom_rate","Vn_fail_rate","Vn_symptom_rate","beta0","beta1","beta_mod_Awp","beta_mod_Aap","beta_mod_An","rho","sigmaSE","lag",init.names)

rproc <- Csnippet("
                  double N, beta, foi, dw, nbirths, Sbirths, Vwp_births, Vap_births;
                  double rate[22], trans[22];
                  
                  // pop size
                  N = S+I+Vwp+Vap+Vn+Awp+Aap+An;
                  
                  // force of infection
                  beta = beta0 * (1 + beta1 * sin(2 * 3.14159 * (t-lag*365.25)/365));
                  foi = (beta*I+beta*beta_mod_Awp*Awp+beta*beta_mod_Aap*Aap+beta*beta_mod_An*An)/N;
                  
                  // white noise (extra-demographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  // S
                  rate[0] = foi*dw/dt;  //infection rate (stochastic)
                  rate[1] = death_rate/365.25; // natural mortality
                  
                  // I
                  rate[2] = rec_rate;  // I recovery rate
                  rate[3] = death_rate/365.25; // natural mortality
                  
                  //Vwp
                  rate[4] = Vwp_fail_rate*Vwp_symptom_rate*foi*dw/dt; // rate at which Vwp becomes I
                  rate[5] = Vwp_fail_rate*(1.0-Vwp_symptom_rate)*foi*dw/dt; // rate at which Vwp becomes Awp
                  rate[6] = Vwp_wane_rate;         // Vwp waning rate
                  rate[7] = death_rate/365.25; // natural mortality
                  
                  //Vap
                  rate[8] = Vap_fail_rate*Vap_symptom_rate*foi*dw/dt; // rate at which Vap becomes I
                  rate[9] = Vap_fail_rate*(1.0-Vap_symptom_rate)*foi*dw/dt; // rate at which Vap becomes Aap
                  rate[10] = Vap_wane_rate;         // Vap waning rate
                  rate[11] = death_rate/365.25; // natural mortality
                  
                  //Vn
                  rate[12] = Vn_fail_rate*Vn_symptom_rate*foi*dw/dt; // rate at which Vn becomes I
                  rate[13] = Vn_fail_rate*(1.0-Vn_symptom_rate)*foi*dw/dt; // rate at which Vn becomes An
                  rate[14] = Vn_wane_rate;         // Vn waning rate
                  rate[15] = death_rate/365.25; // natural mortality
                  
                  //Awp
                  rate[16] = rec_rate;   // Awp recovery rate
                  rate[17] = death_rate/365.25; // natural mortality
                  
                  //Aap
                  rate[18] = rec_rate;   // Aap recovery rate
                  rate[19] = death_rate/365.25; // natural mortality
                  
                  //An
                  rate[20] = rec_rate;   // An recovery rate
                  rate[21] = death_rate/365.25; // natural mortality
                  
                  // Poisson births
                  nbirths = rpois((birth_rate/365.25)*N*dt);
                  Sbirths = nearbyint((1.0-vacc_rate)*nbirths);
                  
                  // Vaccinations
                  if (t<9862) {Vwp_births = nearbyint(vacc_rate*nbirths); Vap_births=0;}
                  if (t>=9862) {Vap_births = nearbyint(vacc_rate*nbirths); Vwp_births=0;} 
                  
                  // transitions between classes
                  reulermultinom(2, S, &rate[0], dt, &trans[0]);
                  reulermultinom(2, I, &rate[2], dt, &trans[2]);
                  reulermultinom(4, Vwp, &rate[4], dt, &trans[4]);
                  reulermultinom(4, Vap, &rate[8], dt, &trans[8]);
                  reulermultinom(4, Vn, &rate[12], dt, &trans[12]);
                  reulermultinom(2, Awp, &rate[16], dt, &trans[16]);
                  reulermultinom(2, Aap, &rate[18], dt, &trans[18]);
                  reulermultinom(2, An, &rate[20], dt, &trans[20]);
                  
                  S += Sbirths + trans[6] + trans[10] + trans[14]- trans[0] - trans[1];
                  I += trans[0] + trans[4] + trans[8] + trans[12]- trans[2] - trans[3];
                  Vwp += Vwp_births - trans[4] - trans[5] - trans[6] - trans[7];
                  Vap += Vap_births - trans[8] - trans[9] - trans[10] - trans[11];
                  Vn += trans[2] + trans[16] + trans[18] + trans [20] - trans[12] - trans[13] - trans[14] - trans[15];
                  Awp += trans[5] - trans[16] - trans[17];
                  Aap += trans[9] - trans[18] - trans[19];
                  An += trans[13] - trans[20] - trans[21];
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[0] + trans[4] + trans[8] + trans[12]; // true incidence
                  
                  if(I<0) {I=0;}
                  ")

rinit <- Csnippet("
                  double m = pop/(S_0+I_0+Vwp_0+Vap_0+Vn_0+Awp_0+Aap_0+An_0);
                  S = nearbyint(m*S_0);
                  I = nearbyint(m*I_0);
                  Vwp = nearbyint(m*Vwp_0);
                  Vap = nearbyint(m*Vap_0);
                  Vn = nearbyint(m*Vn_0);
                  Awp = nearbyint(m*Awp_0);
                  Aap = nearbyint(m*Aap_0);
                  An = nearbyint(m*An_0);
                  W = 0.0;
                  C = 0.0;
                  ")

dmeas <- Csnippet("
                  lik = dpois(incidence,rho*C,give_log);
                  //if (!R_FINITE(lik)) {Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",S,I,Vwp,Vap,Vn,Awp,Aap,An,beta0,beta1,lag,lik,C);}
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
  log=c("Vwp_wane_rate","Vap_wane_rate","Vn_wane_rate","rec_rate","beta0","sigmaSE"),
  logit=c("Vwp_fail_rate","Vwp_symptom_rate","Vap_fail_rate","Vap_symptom_rate","Vn_fail_rate","Vn_symptom_rate","beta_mod_Awp","beta_mod_Aap","beta_mod_An","rho","lag","beta1"),
  barycentric=c("S_0","I_0","Vwp_0","Vap_0","Vn_0","Awp_0","Aap_0","An_0")
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

est.pars<-c("Vwp_wane_rate","Vap_wane_rate","Vn_wane_rate","rec_rate","beta0","sigmaSE",
            "Vwp_fail_rate","Vwp_symptom_rate","Vap_fail_rate","Vap_symptom_rate","Vn_fail_rate","Vn_symptom_rate",
            "beta_mod_Awp","beta_mod_Aap","beta_mod_An","rho","lag","beta1","S_0","I_0","Vwp_0","Vap_0","Vn_0","Awp_0","Aap_0","An_0")

rw.sd=rw.sd(Vwp_wane_rate=0.02,Vap_wane_rate=ifelse(time >= 9862,0.02,0),Vn_wane_rate=0.02,rec_rate=0.02,beta0=0.02,
            sigmaSE=0.02,Vwp_fail_rate=0.02,Vwp_symptom_rate=0.02,Vap_fail_rate=ifelse(time >= 9862,0.02,0),Vap_symptom_rate=ifelse(time >= 9862,0.02,0),
            Vn_fail_rate=0.02,Vn_symptom_rate=0.02,beta_mod_Awp=0.02,beta_mod_Aap=ifelse(time >= 9862,0.02,0),beta_mod_An=0.02,
            rho=0.02,lag=0.02,beta1=ivp(0.02),S_0=ivp(0.02),I_0=ivp(0.02),Vwp_0=ivp(0.02),Vap_0=ivp(0.02),Vn_0=ivp(0.02),Awp_0=ivp(0.02),Aap_0=ivp(0.02),An_0=ivp(0.02))

init.values.mat<-data.frame(S_0=LHS[,"S_0"],I_0=LHS[,"I_0"],Vwp_0=LHS[,"Vwp_0"],Vap_0=LHS[,"Vap_0"],Vn_0=LHS[,"Vn_0"],Awp_0=LHS[,"Awp_0"],Aap_0=LHS[,"Aap_0"],An_0=LHS[,"An_0"])

params.mat<-cbind(data.frame(Vwp_wane_rate = LHS[,"Vwp_wane_rate"],Vap_wane_rate = LHS[,"Vap_wane_rate"], Vn_wane_rate = LHS[,"Vn_wane_rate"], rec_rate = LHS[,"rec_rate"],
          Vwp_fail_rate=LHS[,"Vwp_fail_rate"], Vwp_symptom_rate=LHS[,"Vwp_symptom_rate"], Vap_fail_rate=LHS[,"Vap_fail_rate"], Vap_symptom_rate=LHS[,"Vap_symptom_rate"], 
          Vn_fail_rate=LHS[,"Vn_fail_rate"],Vn_symptom_rate=LHS[,"Vn_symptom_rate"],beta0 = LHS[,"beta0"], beta1 = LHS[,"beta1"],
          beta_mod_Awp = LHS[,"beta_mod_Awp"], beta_mod_Aap = LHS[,"beta_mod_Aap"], beta_mod_An = LHS[,"beta_mod_An"], rho=LHS[,"rho"],
          sigmaSE=LHS[,"sigmaSE"],lag=LHS[,"lag"]),init.values.mat)

