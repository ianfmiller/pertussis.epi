#### build pomp model--test model ###

library(lubridate)
library(pomp)
library(doParallel)

statenames <- c("S","I","A","Vv","Vn","W","C")
t0 <- min(dates)

init.names<-c("S_0","I_0","A_0","Vv_0","Vn_0")
param.names<-c("Vv_wane_rate","Vn_wane_rate","rec_rate","V_fail_rate","V_symptom_rate","beta0","beta1","beta_mod_A","rho","sigmaSE","lag","TDaP_waning_boost",init.names)

rproc <- Csnippet("
                  double N, beta, foi, dw, nbirths, Sbirths, vacc_births, wane_mod_TDaP;
                  double rate[14], trans[14];
                  
                  // pop size
                  N = S+I+A+Vv+Vn;
                  
                  // force of infection
                  beta = beta0 * (1 + beta1 * sin(2 * 3.14159 * (t-lag*365.25)/365));
                  foi = (beta*I+beta*beta_mod_A*A)/N;
                  
                  // white noise (extra-demographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  // TDaP waning scalar
                  if (t<13149) {wane_mod_TDaP=1;}
                  if (t>=13149) {wane_mod_TDaP=TDaP_waning_boost;}
                  
                  //S
                  rate[0] = foi*dw/dt;  //infection rate (stochastic)
                  rate[1] = death_rate/365.25; // natural mortality
                  
                  //I
                  rate[2] = rec_rate;  //recovery rate
                  rate[3] = death_rate/365.25; // natural mortality
                  
                  //A
                  rate[4] = rec_rate;   // recovery rate
                  rate[5] = death_rate/365.25; // natural mortality
                  
                  //Vv
                  rate[6] = V_fail_rate*V_symptom_rate*foi*dw/dt; // rate at which Vv becomes I
                  rate[7] = V_fail_rate*(1.0-V_symptom_rate)*foi*dw/dt; // rate at which Vv becomes A
                  rate[8] = Vv_wane_rate*wane_mod_TDaP; // Vv waning rate
                  rate[9] = death_rate/365.25; // natural mortality
                  
                  //Vn
                  rate[10] = V_fail_rate*V_symptom_rate*foi*dw/dt; // rate at which Vn becomes I
                  rate[11] = V_fail_rate*(1.0-V_symptom_rate)*foi*dw/dt; // rate at which Vn becomes A
                  rate[12] = Vn_wane_rate*wane_mod_TDaP;   // Vn waning rate
                  rate[13] = death_rate/365.25; // natural mortality
                  
                  // Poisson births
                  nbirths = rpois((birth_rate/365.25)*N*dt);
                  Sbirths = nearbyint((1.0-vacc_rate)*nbirths);
                  vacc_births = nearbyint(vacc_rate*nbirths);
                  
                  // transitions between classes
                  reulermultinom(2, S, &rate[0], dt, &trans[0]);
                  reulermultinom(2, I, &rate[2], dt, &trans[2]);
                  reulermultinom(2, A, &rate[4], dt, &trans[4]);
                  reulermultinom(4, Vv, &rate[6], dt, &trans[6]);
                  reulermultinom(4, Vn, &rate[10], dt, &trans[10]);
      
                  
                  S += Sbirths + trans[8] + trans[12] - trans[0] - trans[1];
                  I += trans[0] + trans[6] + trans[10]- trans[2] - trans[3];
                  A += trans[7] + trans[11] - trans[4] - trans[5];
                  Vv += vacc_births - trans[6] - trans[7] - trans[8] - trans[9];
                  Vn += trans[2] + trans[4] - trans[10] - trans[11] - trans[12] - trans[13];
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[0] + trans[6] + trans[10];   // true incidence
                  
                  if(I<0) {I=0;}
                  ")

rinit <- Csnippet("
                  double m = pop/(S_0+I_0+A_0+Vv_0+Vn_0);
                  S = nearbyint(m*S_0);
                  I = nearbyint(m*I_0);
                  A = nearbyint(m*A_0);
                  Vv = nearbyint(m*Vv_0);
                  Vn = nearbyint(m*Vn_0);
                  W = 0.0;
                  C = 0.0;
                  ")

dmeas <- Csnippet("
                  lik = dpois(incidence,rho*C,give_log);
                  //if (!R_FINITE(lik)) {Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",S,I,A,Vv,Vn,beta0,beta1,lag,lik,C);}
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
  log=c("TDaP_waning_boost","Vv_wane_rate","Vn_wane_rate","rec_rate","beta0","sigmaSE"),
  logit=c("V_fail_rate","V_symptom_rate","beta_mod_A","rho","lag","beta1"),
  barycentric=c("S_0","I_0","A_0","Vv_0","Vn_0")
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

est.pars<-c("TDaP_waning_boost","Vv_wane_rate","Vn_wane_rate","rec_rate","beta0","sigmaSE",
            "V_fail_rate","V_symptom_rate","beta_mod_A","rho",
            "lag","beta1","S_0","I_0","A_0","Vv_0","Vn_0")

rw.sd=rw.sd(TDaP_waning_boost=0.02,Vv_wane_rate=0.02,Vn_wane_rate=0.02,rec_rate=0.02,beta0=0.02,sigmaSE=0.02,
            V_fail_rate=0.02,V_symptom_rate=0.02,beta_mod_A=0.02,rho=0.02,lag=0.02,
            beta1=0.02,S_0=0.02,I_0=0.02,A_0=0.02,Vv_0=0.02,Vn_0=0.02)

init.values.mat<-data.frame(S_0=LHS[,"S_0"],I_0=LHS[,"I_0"],A_0=LHS[,"A_0"],Vv_0=LHS[,"Vv_0"],Vn_0=LHS[,"Vn_0"])

params.mat<-cbind(data.frame(TDaP_waning_boost = LHS[,"TDaP_waning_boost"],Vv_wane_rate = LHS[,"Vv_wane_rate"],Vn_wane_rate = LHS[,"Vn_wane_rate"],rec_rate = LHS[,"rec_rate"],
                             V_fail_rate=LHS[,"V_fail_rate"], V_symptom_rate=LHS[,"V_symptom_rate"],
                             beta0 = LHS[,"beta0"], beta1 = LHS[,"beta1"],
                             beta_mod_A = LHS[,"beta_mod_A"], rho=LHS[,"rho"],
                             sigmaSE=LHS[,"sigmaSE"],lag=LHS[,"lag"]),init.values.mat)

