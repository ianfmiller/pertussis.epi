#### build pomp model--test model ###

library(lubridate)
library(pomp)
library(doParallel)

statenames <- c("S","I","Vv","Vn","Av","An","W","C")
t0 <- min(dates)

init.names<-c("S_0","I_0","Vv_0","Vn_0","Av_0","An_0")
param.names<-c("Vv_wane_rate","Vn_wane_rate","rec_rate","Vv_fail_rate","Vv_symptom_rate","Vn_fail_rate","Vn_symptom_rate","beta0","beta1","beta_mod_Av","beta_mod_An","rho","sigmaSE","lag",init.names)

rproc <- Csnippet("
                  double N, beta, foi, dw, nbirths, Sbirths, vacc_births;
                  double rate[16], trans[16];
                  
                  // pop size
                  N = S+I+Vv+Vn+Av+An;
                  
                  // force of infection
                  beta = beta0 * (1 + beta1 * sin(2 * 3.14159 * (t-lag*365.25)/365));
                  foi = (beta*I+beta*beta_mod_Av*Av+beta*beta_mod_An*An)/N;
                  
                  // white noise (extra-demographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  // S
                  rate[0] = foi*dw/dt;  //infection rate (stochastic)
                  rate[1] = death_rate/365.25; // natural mortality
                  
                  // I
                  rate[2] = rec_rate;  // I recovery rate
                  rate[3] = death_rate/365.25; // natural mortality
                  
                  //Vv
                  rate[4] = Vv_fail_rate*Vv_symptom_rate*foi*dw/dt; // rate at which Vv becomes I
                  rate[5] = Vv_fail_rate*(1.0-Vv_symptom_rate)*foi*dw/dt; // rate at which Vv becomes Av
                  rate[6] = Vv_wane_rate;         // Vv waning rate
                  rate[7] = death_rate/365.25; // natural mortality
                  
                  //Vn
                  rate[8] = Vn_fail_rate*Vn_symptom_rate*foi*dw/dt; // rate at which Vn becomes I
                  rate[9] = Vn_fail_rate*(1.0-Vn_symptom_rate)*foi*dw/dt; // rate at which Vn becomes An
                  rate[10] = Vn_wane_rate;         // Vn waning rate
                  rate[11] = death_rate/365.25; // natural mortality
                  
                  //Av
                  rate[12] = rec_rate;   // Awp recovery rate
                  rate[13] = death_rate/365.25; // natural mortality
                  
                  //An
                  rate[14] = rec_rate;   // An recovery rate
                  rate[15] = death_rate/365.25; // natural mortality
                  
                  // Poisson births
                  nbirths = rpois((birth_rate/365.25)*N*dt);
                  Sbirths = nearbyint((1.0-vacc_rate)*nbirths);
                  vacc_births = nearbyint(vacc_rate*nbirths);
                  
                  // transitions between classes
                  reulermultinom(2, S, &rate[0], dt, &trans[0]);
                  reulermultinom(2, I, &rate[2], dt, &trans[2]);
                  reulermultinom(4, Vv, &rate[4], dt, &trans[4]);
                  reulermultinom(4, Vn, &rate[8], dt, &trans[8]);
                  reulermultinom(2, Av, &rate[12], dt, &trans[12]);
                  reulermultinom(2, An, &rate[14], dt, &trans[14]);
                  
                  S += Sbirths + trans[6] + trans[10] - trans[0] - trans[1];
                  I += trans[0] + trans[4] + trans[8] - trans[2] - trans[3];
                  Vv += vacc_births - trans[4] - trans[5] - trans[6] - trans[7];
                  Vn += trans[2] + trans[12] + trans[14]  - trans[8] - trans[9] - trans[10] - trans[11];
                  Av += trans[5] - trans[12] - trans[13];
                  An += trans[9] - trans[14] - trans[15];
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[0] + trans[4] + trans[8]; // true incidence
                  
                  if(I<0) {I=0;}
                  ")

rinit <- Csnippet("
                  double m = pop/(S_0+I_0+Vv_0+Vn_0+Av_0+An_0);
                  S = nearbyint(m*S_0);
                  I = nearbyint(m*I_0);
                  Vv = nearbyint(m*Vv_0);
                  Vn = nearbyint(m*Vn_0);
                  Av = nearbyint(m*Av_0);
                  An = nearbyint(m*An_0);
                  W = 0.0;
                  C = 0.0;
                  ")

dmeas <- Csnippet("
                  lik = dpois(incidence,rho*C,give_log);
                  //if (!R_FINITE(lik)) {Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",S,I,Vv,Vn,Av,An,beta0,beta1,lag,lik,C);}
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
  log=c("Vv_wane_rate","Vn_wane_rate","rec_rate","beta0","sigmaSE"),
  logit=c("Vv_fail_rate","Vv_symptom_rate","Vn_fail_rate","Vn_symptom_rate","beta_mod_Av","beta_mod_An","rho","lag","beta1"),
  barycentric=c("S_0","I_0","Vv_0","Vn_0","Av_0","An_0")
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

est.pars<-c("Vv_wane_rate","Vn_wane_rate","rec_rate","beta0","sigmaSE",
            "Vv_fail_rate","Vv_symptom_rate","Vn_fail_rate","Vn_symptom_rate",
            "beta_mod_Av","beta_mod_An","rho","lag","beta1","S_0","I_0","Vv_0","Vn_0","Av_0","An_0")

rw.sd=rw.sd(Vv_wane_rate=0.02,Vn_wane_rate=0.02,rec_rate=0.02,beta0=0.02,
            sigmaSE=0.02,Vv_fail_rate=0.02,Vv_symptom_rate=0.02,
            Vn_fail_rate=0.02,Vn_symptom_rate=0.02,beta_mod_Av=0.02,beta_mod_An=0.02,
            rho=0.02,lag=0.02,beta1=0.02,S_0=0.02,I_0=0.02,Vv_0=0.02,Vn_0=0.02,Av_0=0.02,An_0=0.02)

init.values.mat<-data.frame(S_0=LHS[,"S_0"],I_0=LHS[,"I_0"],Vv_0=LHS[,"Vv_0"],Vn_0=LHS[,"Vn_0"],Av_0=LHS[,"Av_0"],An_0=LHS[,"An_0"])

params.mat<-cbind(data.frame(Vv_wane_rate = LHS[,"Vv_wane_rate"], Vn_wane_rate = LHS[,"Vn_wane_rate"], rec_rate = LHS[,"rec_rate"],
                             Vv_fail_rate=LHS[,"Vv_fail_rate"], Vv_symptom_rate=LHS[,"Vv_symptom_rate"],
                             Vn_fail_rate=LHS[,"Vn_fail_rate"],Vn_symptom_rate=LHS[,"Vn_symptom_rate"],beta0 = LHS[,"beta0"], beta1 = LHS[,"beta1"],
                             beta_mod_Av = LHS[,"beta_mod_Av"],  beta_mod_An = LHS[,"beta_mod_An"], rho=LHS[,"rho"],
                             sigmaSE=LHS[,"sigmaSE"],lag=LHS[,"lag"]),init.values.mat)

