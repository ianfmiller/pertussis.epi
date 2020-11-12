#### build pomp model--test model ###

library(lubridate)
library(pomp)
library(doParallel)

statenames <- c("S","I","Vnwp","Vap","Anwp","Aap","W","C")
t0 <- min(dates)

init.names<-c("S_0","I_0","Vnwp_0","Vap_0","Anwp_0","Aap_0")
param.names<-c("Vnwp_wane_rate","Vap_wane_rate","rec_rate","Vnwp_fail_rate","Vnwp_symptom_rate","Vap_fail_rate","Vap_symptom_rate","beta0","beta1","beta_mod_Anwp","beta_mod_Aap","rho","sigmaSE","lag",init.names)

rproc <- Csnippet("
                  double N, beta, foi, dw, nbirths, Sbirths, Vwp_births, Vap_births;
                  double rate[16], trans[16];
                  
                  // pop size
                  N = S+I+Vnwp+Vap+Anwp+Aap;
                  
                  // force of infection
                  beta = beta0 * (1 + beta1 * sin(2 * 3.14159 * (t-lag*365.25)/365));
                  foi = (beta*I+beta*beta_mod_Anwp*Anwp+beta*beta_mod_Aap*Aap)/N;
                  
                  // white noise (extra-demographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  // S
                  rate[0] = foi*dw/dt;  //infection rate (stochastic)
                  rate[1] = death_rate/365.25; // natural mortality
                  
                  // I
                  rate[2] = rec_rate;  // I recovery rate
                  rate[3] = death_rate/365.25; // natural mortality
                  
                  //Vnwp
                  rate[4] = Vnwp_fail_rate*Vnwp_symptom_rate*foi*dw/dt; // rate at which Vnwp becomes I
                  rate[5] = Vnwp_fail_rate*(1.0-Vnwp_symptom_rate)*foi*dw/dt; // rate at which Vnwp becomes Anwp
                  rate[6] = Vnwp_wane_rate;         // Vnwp waning rate
                  rate[7] = death_rate/365.25; // natural mortality
                  
                  //Vap
                  rate[8] = Vap_fail_rate*Vap_symptom_rate*foi*dw/dt; // rate at which Vap becomes I
                  rate[9] = Vap_fail_rate*(1.0-Vap_symptom_rate)*foi*dw/dt; // rate at which Vap becomes Aap
                  rate[10] = Vap_wane_rate;         // Vap waning rate
                  rate[11] = death_rate/365.25; // natural mortality
                  
                  //Anwp
                  rate[12] = rec_rate;   // Anwp recovery rate
                  rate[13] = death_rate/365.25; // natural mortality
                  
                  //Aap
                  rate[14] = rec_rate;   // Aap recovery rate
                  rate[15] = death_rate/365.25; // natural mortality
                  
                  // Poisson births
                  nbirths = rpois((birth_rate/365.25)*N*dt);
                  Sbirths = nearbyint((1.0-vacc_rate)*nbirths);
                  
                  // Vaccinations
                  if (t<9862) 
                  {
                  Vwp_births = nearbyint(vacc_rate*nbirths); 
                  Vap_births=0;
                  }
                  
                  if (t>=9862) 
                  {
                  Vap_births = nearbyint(vacc_rate*nbirths);
                  Vwp_births=0;
                  } 
                  
                  // transitions between classes
                  reulermultinom(2, S, &rate[0], dt, &trans[0]);
                  reulermultinom(2, I, &rate[2], dt, &trans[2]);
                  reulermultinom(4, Vnwp, &rate[4], dt, &trans[4]);
                  reulermultinom(4, Vap, &rate[8], dt, &trans[8]);
                  reulermultinom(2, Anwp, &rate[12], dt, &trans[12]);
                  reulermultinom(2, Aap, &rate[14], dt, &trans[14]);

                  S += Sbirths + trans[6] + trans[10] - trans[0] - trans[1];
                  I += trans[0] + trans[4] + trans[8] - trans[2] - trans[3];
                  Vnwp += Vwp_births + trans[2] + trans[12] + trans[14] - trans[4] - trans[5] - trans[6] - trans[7];
                  Vap += Vap_births - trans[8] - trans[9] - trans[10] - trans[11];
                  Anwp += trans[5] - trans[12] - trans[13];
                  Aap += trans[9] - trans[14] - trans[15];
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[0] + trans[4] + trans[8]; // true incidence
                  
                  if(I<0) {I=0;}
                  ")

rinit <- Csnippet("
                  double m = pop/(S_0+I_0+Vnwp_0+Vap_0+Anwp_0+Aap_0);
                  S = nearbyint(m*S_0);
                  I = nearbyint(m*I_0);
                  Vnwp = nearbyint(m*Vnwp_0);
                  Vap = nearbyint(m*Vap_0);
                  Anwp = nearbyint(m*Anwp_0);
                  Aap = nearbyint(m*Aap_0);
                  W = 0.0;
                  C = 0.0;
                  ")

dmeas <- Csnippet("
                  lik = dpois(incidence,rho*C,give_log);
                  //if (!R_FINITE(lik)) {Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",S,I,Vnwp,Vap,Anwp,Aap,beta0,beta1,lag,lik,C);}
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
  log=c("Vnwp_wane_rate","Vap_wane_rate","rec_rate","beta0","sigmaSE"),
  logit=c("Vnwp_fail_rate","Vnwp_symptom_rate","Vap_fail_rate","Vap_symptom_rate","beta_mod_Anwp","beta_mod_Aap","rho","lag","beta1"),
  barycentric=c("S_0","I_0","Vnwp_0","Vap_0","Anwp_0","Aap_0")
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

est.pars<-c("Vnwp_wane_rate","Vap_wane_rate","rec_rate","beta0","sigmaSE",
            "Vnwp_fail_rate","Vnwp_symptom_rate","Vap_fail_rate","Vap_symptom_rate",
            "beta_mod_Anwp","beta_mod_Aap","rho","lag","beta1","S_0","I_0","Vnwp_0","Vap_0","Anwp_0","Aap_0")

rw.sd=rw.sd(Vnwp_wane_rate=0.02,Vap_wane_rate=0.02,rec_rate=0.02,beta0=0.02,
            sigmaSE=0.02,Vnwp_fail_rate=0.02,Vnwp_symptom_rate=0.02,Vap_fail_rate=0.02,Vap_symptom_rate=0.02,
            beta_mod_Anwp=0.02,beta_mod_Aap=0.02,
            rho=0.02,lag=0.02,beta1=0.02,S_0=0.02,I_0=0.02,Vnwp_0=0.02,Vap_0=0.02,Anwp_0=0.02,Aap_0=0.02)

init.values.mat<-data.frame(S_0=LHS[,"S_0"],I_0=LHS[,"I_0"],Vnwp_0=LHS[,"Vnwp_0"],Vap_0=LHS[,"Vap_0"],Anwp_0=LHS[,"Anwp_0"],Aap_0=LHS[,"Aap_0"])

params.mat<-cbind(data.frame(Vnwp_wane_rate = LHS[,"Vnwp_wane_rate"],Vap_wane_rate = LHS[,"Vap_wane_rate"], rec_rate = LHS[,"rec_rate"],
                             Vnwp_fail_rate=LHS[,"Vnwp_fail_rate"], Vnwp_symptom_rate=LHS[,"Vnwp_symptom_rate"], Vap_fail_rate=LHS[,"Vap_fail_rate"], Vap_symptom_rate=LHS[,"Vap_symptom_rate"], 
                             beta0 = LHS[,"beta0"], beta1 = LHS[,"beta1"],
                             beta_mod_Anwp = LHS[,"beta_mod_Anwp"], beta_mod_Aap = LHS[,"beta_mod_Aap"], rho=LHS[,"rho"],
                             sigmaSE=LHS[,"sigmaSE"],lag=LHS[,"lag"]),init.values.mat)

