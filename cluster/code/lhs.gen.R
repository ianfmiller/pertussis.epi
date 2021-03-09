### generalte/load LHS ###
setwd(lhs.dir)
LHS.file.name<-paste(loc,model,subset.data,smooth.interval,"LHS.RDS",sep=".")
if (file.exists(LHS.file.name))
{
  LHS<-readRDS(LHS.file.name)
} else
{
  n.points<-n.initial
  
  if (model=="waning.diff.protection.diff.symptoms.diff.transmission.diff")
  {
    lower.param.bounds<-c(
      Vn_wane_rate = 1/(365*1),
      Vwp_wane_rate = 1/(365*1), 
      Vap_wane_rate = 1/(365*1),
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      Vwp_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vap_fail_rate=0,
      Vap_symptom_rate=0, 
      naive_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0.0,
      beta_mod_An = 0.0, 
      beta_mod_Awp = 0.0, 
      beta_mod_Aap = 0.0, 
      beta_mod_A = 0.0,
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      E_0=0.0,
      En_0=0.0,
      Ewp_0=0.0,
      I_0=0.0,
      A_0=0.0,
      An_0=0.0,
      Awp_0=0.0,
      Vn_0=0.0,
      Vwp_0=0.0)
    
    upper.param.bounds<-c(
      Vn_wane_rate = 1/(365*500),
      Vwp_wane_rate = 1/(365*500), 
      Vap_wane_rate = 1/(365*500),
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      Vwp_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vap_fail_rate=1,
      Vap_symptom_rate=1, 
      naive_symptom_rate=1,
      beta0 = 1, 
      beta1 = 1,
      beta_mod_An = 1, 
      beta_mod_Awp = 1, 
      beta_mod_Aap = 1, 
      beta_mod_A = 1,
      rho=1,
      sigmaSE=0.25,
      lag=1,
      S_0=1,
      E_0=.1,
      En_0=.1,
      Ewp_0=.1,
      I_0=.1,
      A_0=.1,
      An_0=.1,
      Awp_0=.1,
      Vn_0=1,
      Vwp_0=1)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="waning.same.protection.diff.symptoms.diff.transmission.diff")
  {
    lower.param.bounds<-c(
      Vn_wane_rate = 1/(365*1),
      V_wane_rate = 1/(365*1), 
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      Vwp_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vap_fail_rate=0,
      Vap_symptom_rate=0, 
      naive_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0.0,
      beta_mod_An = 0.0, 
      beta_mod_Awp = 0.0, 
      beta_mod_Aap = 0.0, 
      beta_mod_A = 0.0,
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      E_0=0.0,
      En_0=0.0,
      Ewp_0=0.0,
      I_0=0.0,
      A_0=0.0,
      An_0=0.0,
      Awp_0=0.0,
      Vn_0=0.0,
      Vwp_0=0.0)
    
    upper.param.bounds<-c(
      Vn_wane_rate = 1/(365*500),
      V_wane_rate = 1/(365*500), 
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      Vwp_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vap_fail_rate=1,
      Vap_symptom_rate=1, 
      naive_symptom_rate=1,
      beta0 = 1, 
      beta1 = 1,
      beta_mod_An = 1, 
      beta_mod_Awp = 1, 
      beta_mod_Aap = 1, 
      beta_mod_A = 1,
      rho=1,
      sigmaSE=0.25,
      lag=1,
      S_0=1,
      E_0=.1,
      En_0=.1,
      Ewp_0=.1,
      I_0=.1,
      A_0=.1,
      An_0=.1,
      Awp_0=.1,
      Vn_0=1,
      Vwp_0=1)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="waning.diff.protection.same.symptoms.diff.transmission.diff")
  {
    lower.param.bounds<-c(
      Vn_wane_rate = 1/(365*1),
      Vwp_wane_rate = 1/(365*1), 
      Vap_wane_rate = 1/(365*1),
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      V_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vap_symptom_rate=0, 
      naive_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0.0,
      beta_mod_An = 0.0, 
      beta_mod_Awp = 0.0, 
      beta_mod_Aap = 0.0, 
      beta_mod_A = 0.0,
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      E_0=0.0,
      En_0=0.0,
      Ewp_0=0.0,
      I_0=0.0,
      A_0=0.0,
      An_0=0.0,
      Awp_0=0.0,
      Vn_0=0.0,
      Vwp_0=0.0)
    
    upper.param.bounds<-c(
      Vn_wane_rate = 1/(365*500),
      Vwp_wane_rate = 1/(365*500), 
      Vap_wane_rate = 1/(365*500),
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      V_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vap_symptom_rate=1, 
      naive_symptom_rate=1,
      beta0 = 1, 
      beta1 = 1,
      beta_mod_An = 1, 
      beta_mod_Awp = 1, 
      beta_mod_Aap = 1, 
      beta_mod_A = 1,
      rho=1,
      sigmaSE=0.25,
      lag=1,
      S_0=1,
      E_0=.1,
      En_0=.1,
      Ewp_0=.1,
      I_0=.1,
      A_0=.1,
      An_0=.1,
      Awp_0=.1,
      Vn_0=1,
      Vwp_0=1)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="waning.diff.protection.diff.symptoms.same.transmission.diff")
  {
    lower.param.bounds<-c(
      Vn_wane_rate = 1/(365*1),
      Vwp_wane_rate = 1/(365*1), 
      Vap_wane_rate = 1/(365*1),
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      Vwp_fail_rate=0, 
      V_symptom_rate=0, 
      Vap_fail_rate=0,
      naive_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0.0,
      beta_mod_An = 0.0, 
      beta_mod_Awp = 0.0, 
      beta_mod_Aap = 0.0, 
      beta_mod_A = 0.0,
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      E_0=0.0,
      En_0=0.0,
      Ewp_0=0.0,
      I_0=0.0,
      A_0=0.0,
      An_0=0.0,
      Awp_0=0.0,
      Vn_0=0.0,
      Vwp_0=0.0)
    
    upper.param.bounds<-c(
      Vn_wane_rate = 1/(365*500),
      Vwp_wane_rate = 1/(365*500), 
      Vap_wane_rate = 1/(365*500),
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      Vwp_fail_rate=1, 
      V_symptom_rate=1, 
      Vap_fail_rate=1,
      naive_symptom_rate=1,
      beta0 = 1, 
      beta1 = 1,
      beta_mod_An = 1, 
      beta_mod_Awp = 1, 
      beta_mod_Aap = 1, 
      beta_mod_A = 1,
      rho=1,
      sigmaSE=0.25,
      lag=1,
      S_0=1,
      E_0=.1,
      En_0=.1,
      Ewp_0=.1,
      I_0=.1,
      A_0=.1,
      An_0=.1,
      Awp_0=.1,
      Vn_0=1,
      Vwp_0=1)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="waning.diff.protection.diff.symptoms.diff.transmission.same")
  {
    lower.param.bounds<-c(
      Vn_wane_rate = 1/(365*1),
      Vwp_wane_rate = 1/(365*1), 
      Vap_wane_rate = 1/(365*1),
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      Vwp_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vap_fail_rate=0,
      Vap_symptom_rate=0, 
      naive_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0.0,
      beta_mod_An = 0.0, 
      beta_mod_Av = 0.0, 
      beta_mod_A = 0.0,
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      E_0=0.0,
      En_0=0.0,
      Ewp_0=0.0,
      I_0=0.0,
      A_0=0.0,
      An_0=0.0,
      Av_0=0.0,
      Vn_0=0.0,
      Vwp_0=0.0)
    
    upper.param.bounds<-c(
      Vn_wane_rate = 1/(365*500),
      Vwp_wane_rate = 1/(365*500), 
      Vap_wane_rate = 1/(365*500),
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      Vwp_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vap_fail_rate=1,
      Vap_symptom_rate=1, 
      naive_symptom_rate=1,
      beta0 = 1, 
      beta1 = 1,
      beta_mod_An = 1, 
      beta_mod_Av = 1, 
      beta_mod_A = 1,
      rho=1,
      sigmaSE=0.25,
      lag=1,
      S_0=1,
      E_0=.1,
      En_0=.1,
      Ewp_0=.1,
      I_0=.1,
      A_0=.1,
      An_0=.1,
      Av_0=.1,
      Vn_0=1,
      Vwp_0=1)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="waning.same.protection.same.symptoms.diff.transmission.diff")
  {
    lower.param.bounds<-c(
      Vn_wane_rate = 1/(365*1),
      V_wane_rate = 1/(365*1), 
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      V_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vap_symptom_rate=0, 
      naive_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0.0,
      beta_mod_An = 0.0, 
      beta_mod_Awp = 0.0, 
      beta_mod_Aap = 0.0, 
      beta_mod_A = 0.0,
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      E_0=0.0,
      En_0=0.0,
      Ewp_0=0.0,
      I_0=0.0,
      A_0=0.0,
      An_0=0.0,
      Awp_0=0.0,
      Vn_0=0.0,
      Vwp_0=0.0)
    
    upper.param.bounds<-c(
      Vn_wane_rate = 1/(365*500),
      V_wane_rate = 1/(365*500), 
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      V_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vap_symptom_rate=1, 
      naive_symptom_rate=1,
      beta0 = 1, 
      beta1 = 1,
      beta_mod_An = 1, 
      beta_mod_Awp = 1, 
      beta_mod_Aap = 1, 
      beta_mod_A = 1,
      rho=1,
      sigmaSE=0.25,
      lag=1,
      S_0=1,
      E_0=.1,
      En_0=.1,
      Ewp_0=.1,
      I_0=.1,
      A_0=.1,
      An_0=.1,
      Awp_0=.1,
      Vn_0=1,
      Vwp_0=1)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="waning.same.protection.diff.symptoms.same.transmission.diff")
  {
    lower.param.bounds<-c(
      Vn_wane_rate = 1/(365*1),
      V_wane_rate = 1/(365*1), 
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      Vwp_fail_rate=0, 
      V_symptom_rate=0, 
      Vap_fail_rate=0,
      naive_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0.0,
      beta_mod_An = 0.0, 
      beta_mod_Awp = 0.0, 
      beta_mod_Aap = 0.0, 
      beta_mod_A = 0.0,
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      E_0=0.0,
      En_0=0.0,
      Ewp_0=0.0,
      I_0=0.0,
      A_0=0.0,
      An_0=0.0,
      Awp_0=0.0,
      Vn_0=0.0,
      Vwp_0=0.0)
    
    upper.param.bounds<-c(
      Vn_wane_rate = 1/(365*500),
      V_wane_rate = 1/(365*500), 
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      Vwp_fail_rate=1, 
      V_symptom_rate=1, 
      Vap_fail_rate=1,
      naive_symptom_rate=1,
      beta0 = 1, 
      beta1 = 1,
      beta_mod_An = 1, 
      beta_mod_Awp = 1, 
      beta_mod_Aap = 1, 
      beta_mod_A = 1,
      rho=1,
      sigmaSE=0.25,
      lag=1,
      S_0=1,
      E_0=.1,
      En_0=.1,
      Ewp_0=.1,
      I_0=.1,
      A_0=.1,
      An_0=.1,
      Awp_0=.1,
      Vn_0=1,
      Vwp_0=1)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="waning.diff.protection.diff.symptoms.diff.transmission.same")
  {
    lower.param.bounds<-c(
      Vn_wane_rate = 1/(365*1),
      V_wane_rate = 1/(365*1), 
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      Vwp_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vap_fail_rate=0,
      Vap_symptom_rate=0, 
      naive_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0.0,
      beta_mod_An = 0.0, 
      beta_mod_Av = 0.0, 
      beta_mod_A = 0.0,
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      E_0=0.0,
      En_0=0.0,
      Ewp_0=0.0,
      I_0=0.0,
      A_0=0.0,
      An_0=0.0,
      Av_0=0.0,
      Vn_0=0.0,
      Vwp_0=0.0)
    
    upper.param.bounds<-c(
      Vn_wane_rate = 1/(365*500),
      V_wane_rate = 1/(365*500), 
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      Vwp_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vap_fail_rate=1,
      Vap_symptom_rate=1, 
      naive_symptom_rate=1,
      beta0 = 1, 
      beta1 = 1,
      beta_mod_An = 1, 
      beta_mod_Av = 1, 
      beta_mod_A = 1,
      rho=1,
      sigmaSE=0.25,
      lag=1,
      S_0=1,
      E_0=.1,
      En_0=.1,
      Ewp_0=.1,
      I_0=.1,
      A_0=.1,
      An_0=.1,
      Av_0=.1,
      Vn_0=1,
      Vwp_0=1)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
}
setwd(base.dir)

  
  