#### build pomp model--test model ###

library(lubridate)
library(pomp)
library(doParallel)

statenames <- c("S","E","En","Ewp","Eap","I","A","An","Awp","Aap","Vn","Vwp","Vap","W","C")
t0 <- min(dates)

init.param.names<-c("S_0","E_0","En_0","Ewp_0","I_0","A_0","An_0","Awp_0","Vn_0","Vwp_0")
param.names<-c("Vn_wane_rate","Vwp_wane_rate","Vap_wane_rate","Vn_fail_rate","Vn_symptom_rate","Vwp_fail_rate","Vwp_symptom_rate","Vap_fail_rate","Vap_symptom_rate","naive_symptom_rate","beta0","beta1","beta_mod_An","beta_mod_Awp","beta_mod_Aap","beta_mod_A","rho","sigmaSE","lag",init.param.names)

rproc <- Csnippet("
                  double N, beta, latent_trans_rate, rec_rate, foi, dw, nbirths, Sbirths, Vwp_births, Vap_births;
                  double rate[33], trans[33];
                  
                  // fixed vars vals from https://stm.sciencemag.org/content/scitransmed/suppl/2018/03/26/10.434.eaaj1748.DC1/aaj1748_SM.pdf
                  
                  latent_trans_rate = 0.125;
                  rec_rate = 0.0666;
                  
                  // pop size
                  N = S+E+En+Ewp+Eap+I+A+An+Awp+Aap+Vn+Vwp+Vap;
                  
                  // force of infection
                  beta = beta0 * (1 + beta1 * sin(2 * 3.14159 * (t-lag*365.25)/365));
                  foi = (beta*I+beta*beta_mod_Awp*Awp+beta*beta_mod_Aap*Aap+beta*beta_mod_An*An+beta*beta_mod_A*A)/N;
                  
                  // white noise (extra-demographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  // S
                  rate[0] = foi*dw/dt;  // rate of transition from S to E (dw/dt is stochastic term)
                  rate[1] = death_rate/365.25; // natural mortality
                  
                  // E
                  rate[2] = (1-naive_symptom_rate)*latent_trans_rate; // rate of transition from E to A
                  rate[3] = naive_symptom_rate*latent_trans_rate; // rate of transition from E to I
                  rate[4] = death_rate/365.25; // natural mortality
                  
                  // En
                  rate[5] = (1-Vn_symptom_rate)*latent_trans_rate; // rate of transition from En to An
                  rate[6] = Vn_symptom_rate*latent_trans_rate; // rate of transition from En to I
                  rate[7] = death_rate/365.25; // natural mortality
                  
                  // Ewp
                  rate[8] = (1-Vwp_symptom_rate)*latent_trans_rate; // rate of transition from Ewp to Awp
                  rate[9] = Vwp_symptom_rate*latent_trans_rate; // rate of transition from Ewp to I
                  rate[10] = death_rate/365.25; // natural mortality
                  
                  // Eap
                  rate[11] = (1-Vap_symptom_rate)*latent_trans_rate; // rate of transition from Ewp to Aap
                  rate[12] = Vap_symptom_rate*latent_trans_rate; // rate of transition from Ewp to I
                  rate[13] = death_rate/365.25; // natural mortality
                  
                  // I
                  rate[14] = rec_rate;  // I recovery rate
                  rate[15] = death_rate/365.25; // natural mortality
                  
                  // A
                  rate[16] = rec_rate;   // A recovery rate
                  rate[17] = death_rate/365.25; // natural mortality
                  
                  // An
                  rate[18] = rec_rate;   // An recovery rate
                  rate[19] = death_rate/365.25; // natural mortality
                  
                  // Awp
                  rate[20] = rec_rate;   // Awp recovery rate
                  rate[21] = death_rate/365.25; // natural mortality
                  
                  // Aap
                  rate[22] = rec_rate;   // Aap recovery rate
                  rate[23] = death_rate/365.25; // natural mortality
                  
                   // Vn
                  rate[24] = Vn_fail_rate*foi*dw/dt; // rate at which Vn becomes En
                  rate[25] = Vn_wane_rate;         // Vn waning rate
                  rate[26] = death_rate/365.25; // natural mortality
                  
                  // Vwp
                  rate[27] = Vwp_fail_rate*foi*dw/dt; // rate at which Vwp becomes Ewp
                  rate[28] = Vwp_wane_rate;         // Vwp waning rate
                  rate[29] = death_rate/365.25; // natural mortality
                  
                  // Vap
                  rate[30] = Vap_fail_rate*foi*dw/dt; // rate at which Vap becomes Eap
                  rate[31] = Vap_wane_rate;         // Vap waning rate
                  rate[32] = death_rate/365.25; // natural mortality
                  
                  // Poisson births
                  nbirths = rpois((birth_rate/365.25)*N*dt);
                  Sbirths = nearbyint((1.0-vacc_rate)*nbirths);
                  
                  // Vaccinations
                  if (t<9862) {Vwp_births = nearbyint(vacc_rate*nbirths); Vap_births=0;}
                  if (t>=9862) {Vap_births = nearbyint(vacc_rate*nbirths); Vwp_births=0;} 
                  
                  // transitions between classes
                  reulermultinom(2, S, &rate[0], dt, &trans[0]);
                  reulermultinom(3, E, &rate[2], dt, &trans[2]);
                  reulermultinom(3, En, &rate[5], dt, &trans[5]);
                  reulermultinom(3, Ewp, &rate[8], dt, &trans[8]);
                  reulermultinom(3, Eap, &rate[11], dt, &trans[11]);
                  reulermultinom(2, I, &rate[14], dt, &trans[14]);
                  reulermultinom(2, A, &rate[16], dt, &trans[16]);
                  reulermultinom(2, An, &rate[18], dt, &trans[18]);
                  reulermultinom(2, Awp, &rate[20], dt, &trans[20]);
                  reulermultinom(2, Aap, &rate[22], dt, &trans[22]);
                  reulermultinom(3, Vn, &rate[24], dt, &trans[24]);
                  reulermultinom(3, Vwp, &rate[27], dt, &trans[27]);
                  reulermultinom(3, Vap, &rate[30], dt, &trans[30]);
                  
                  S += Sbirths + trans[25] + trans[28] + trans[31]- trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3] - trans[4];
                  En += trans[24] - trans[5] - trans[6] - trans[7];
                  Ewp += trans[27] - trans[8] - trans[9] - trans[10];
                  Eap += trans[30] - trans[11] - trans[12] - trans[13];
                  I += trans[3] + trans[6] + trans[9] + trans[12]- trans[14] - trans[15];
                  A += trans[2] - trans[16] - trans[17];
                  An += trans[5] - trans[18] - trans[19];
                  Awp += trans[8] - trans[20] - trans[21];
                  Aap += trans[11] - trans[22] - trans[23];
                  Vn += trans[14] + trans[16] + trans[18] + trans [20] + trans[22] - trans[24] - trans[25] - trans[26];
                  Vwp += Vwp_births - trans[27] - trans[28] - trans[29];
                  Vap += Vap_births - trans[30] - trans[31] - trans[32];

                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[3] + trans[6] + trans[9] + trans[12]; // true incidence
                  
                  if(I<0) {I=0;}
                  ")

rinit <- Csnippet("
                  double m = pop/(S_0+E_0+En_0+Ewp_0+I_0+A_0+An_0+Awp_0+Vn_0+Vwp_0);
                  S = nearbyint(m*S_0);
                  E = nearbyint(m*E_0);
                  En = nearbyint(En_0);
                  Ewp = nearbyint(m*Ewp_0);
                  Eap = 0.0;
                  I = nearbyint(m*I_0);
                  A = nearbyint(m*A_0);
                  An = nearbyint(m*An_0);
                  Awp = nearbyint(m*Awp_0);
                  Aap = 0.0;
                  Vn = nearbyint(m*Vn_0);
                  Vwp = nearbyint(m*Vwp_0);
                  Vap = 0.0;
                  W = 0.0;
                  C = 0.0;
                  ")

dmeas <- Csnippet("
                  lik = dpois(incidence,rho*C,give_log);
                  //if (!R_FINITE(lik)) {Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",S,E,En,Ewp,Eap,I,A,An,Awp,Aap,Vn,Vwp,Vap,beta0,beta1,lag,lik,C);}
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
  log=c("Vn_wane_rate","Vwp_wane_rate","Vap_wane_rate","beta0","sigmaSE"),
  logit=c("naive_symptom_rate","Vn_fail_rate","Vn_symptom_rate","Vwp_fail_rate","Vwp_symptom_rate","Vap_fail_rate","Vap_symptom_rate","beta_mod_A","beta_mod_An","beta_mod_Awp","beta_mod_Aap","rho","lag","beta1"),
  barycentric=c("S_0","E_0","En_0","Ewp_0","I_0","Vn_0","Vwp_0","A","An_0","Awp_0")
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

est.pars<-c("Vn_wane_rate","Vwp_wane_rate","Vap_wane_rate","Vn_fail_rate","Vn_symptom_rate","Vwp_fail_rate","Vwp_symptom_rate","Vap_fail_rate","Vap_symptom_rate","naive_symptom_rate",
            "beta0","beta1","beta_mod_An","beta_mod_Awp","beta_mod_Aap","beta_mod_A","rho","sigmaSE","lag",
            "S_0","E_0","En_0","Ewp_0","I_0","Vn_0","Vwp_0","A","An_0","Awp_0")

rw.sd=rw.sd(Vn_wane_rate=0.01,Vwp_wane_rate=0.01,Vap_wane_rate=ifelse(time >= 9862,0.01,0),Vn_fail_rate=0.01,Vn_symptom_rate=0.01,Vwp_fail_rate=0.01,Vwp_symptom_rate=0.01,Vap_fail_rate=ifelse(time >= 9862,0.01,0),Vap_symptom_rate=ifelse(time >= 9862,0.01,0),naive_symptom_rate=0.01,
            beta0=0.01,beta1=0.01,beta_mod_An=0.01,beta_mod_Awp=0.01,beta_mod_Aap=ifelse(time >= 9862,0.01,0),beta_mod_A=0.01,rho=0.01,sigmaSE=0.01,lag=0.01,
            S_0=ivp(0.01),E_0=ivp(0.01),En_0=ivp(0.01),Ewp_0=ivp(0.01),I_0=ivp(0.01),A_0=ivp(0.01),An_0=ivp(0.01),Awp_0=ivp(0.01),Vn_0=ivp(0.01),Vwp_0=ivp(0.01))

init.values.mat<-data.frame(S_0=LHS[,"S_0"],
                            E_0=LHS[,"E_0"],En_0=LHS[,"En_0"],Ewp_0=LHS[,"Ewp_0"],
                            I_0=LHS[,"I_0"],A_0=LHS[,"A_0"],An_0=LHS[,"An_0"],Awp_0=LHS[,"Awp_0"],
                            Vn_0=LHS[,"Vn_0"],Vwp_0=LHS[,"Vwp_0"])

params.mat<-cbind(
  data.frame(
    Vn_wane_rate=LHS[,"Vn_wane_rate"],Vwp_wane_rate = LHS[,"Vwp_wane_rate"],Vap_wane_rate = LHS[,"Vap_wane_rate"], Vn_fail_rate=LHS[,"Vn_fail_rate"],Vn_symptom_rate=LHS[,"Vn_symptom_rate"],
    Vwp_fail_rate=LHS[,"Vwp_fail_rate"],Vwp_symptom_rate=LHS[,"Vwp_symptom_rate"],Vap_fail_rate=LHS[,"Vap_fail_rate"], Vap_symptom_rate=LHS[,"Vap_symptom_rate"],naive_symptom_rate=LHS[,"naive_symptom_rate"],
    beta0 = LHS[,"beta0"], beta1 = LHS[,"beta1"],beta_mod_An = LHS[,"beta_mod_An"],beta_mod_Awp = LHS[,"beta_mod_Awp"], beta_mod_Aap = LHS[,"beta_mod_Aap"],beta_mod_A=LHS[,"beta_mod_A"],rho=LHS[,"rho"],
    sigmaSE=LHS[,"sigmaSE"],lag=LHS[,"lag"]
    ),
  init.values.mat)

