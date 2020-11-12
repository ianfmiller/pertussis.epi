### generalte/load LHS ###
setwd(lhs.dir)
LHS.file.name<-paste(loc,model,subset.data,smooth.interval,"LHS.RDS",sep=".")
if (file.exists(LHS.file.name))
{
  LHS<-readRDS(LHS.file.name)
} else
{
  n.points<-n.initial
  if (model=="test.stoch")
  {
    lower.param.bounds<-c(
      Vwp_wane_rate = 1/(365*1), 
      Vn_wane_rate = 1/(365*1), 
      rec_rate = 1/5,
      Vwp_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_Awp = 0.0, 
      beta_mod_An = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      Vwp_0=0.0,
      Vn_0=0.0,
      Awp_0=0.0,
      An_0=0.0)
    
    upper.param.bounds<-c(
      Vwp_wane_rate = 1/(365*50), 
      Vn_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      Vwp_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_Awp = 1, 
      beta_mod_An = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      Vwp_0=1.0,
      Vn_0=1.0,
      Awp_0=1.0,
      An_0=1.0)
    
      LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
      saveRDS(LHS,file=LHS.file.name)
  }
  if (model=="all.equal")
  {
    lower.param.bounds<-c(
      V_wane_rate = 1/(365*1), 
      rec_rate = 1/5,
      V_fail_rate=0, 
      V_symptom_rate=0, 
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_A = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      A_0=0.0,
      V_0=0.0)
    
    upper.param.bounds<-c(
      V_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      V_fail_rate=1, 
      V_symptom_rate=1, 
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_A = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      A_0=1.0,
      V_0=1.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  if (model=="all.equal.booster")
  {
    lower.param.bounds<-c(
      Vv_wane_rate = 1/(365*1), 
      Vn_wane_rate = 1/(365*1), 
      rec_rate = 1/5,
      V_fail_rate=0, 
      V_symptom_rate=0, 
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_A = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      A_0=0.0,
      Vv_0=0.0,
      Vn_0=0.0,
      TDaP_waning_boost=0.0)
    
    upper.param.bounds<-c(
      Vv_wane_rate = 1/(365*50), 
      Vn_wane_rate = 1/(365*50),
      rec_rate = 1/31,
      V_fail_rate=1, 
      V_symptom_rate=1, 
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_A = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      A_0=1.0,
      Vv_0=1.0,
      Vn_0=1.0,
      TDaP_waning_boost=5.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="Vwp.equal.Vap")
  {
    lower.param.bounds<-c(
      Vv_wane_rate = 1/(365*1), 
      Vn_wane_rate = 1/(365*1), 
      rec_rate = 1/5,
      Vv_fail_rate=0, 
      Vv_symptom_rate=0, 
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_Av = 0.0, 
      beta_mod_An = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      Vv_0=0.0,
      Vn_0=0.0,
      Av_0=0.0,
      An_0=0.0)
    
    upper.param.bounds<-c(
      Vv_wane_rate = 1/(365*50), 
      Vn_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      Vv_fail_rate=1, 
      Vv_symptom_rate=1, 
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_Av = 1, 
      beta_mod_An = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      Vv_0=1.0,
      Vn_0=1.0,
      Av_0=1.0,
      An_0=1.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="Vwp.equal.Vap.booster")
  {
    lower.param.bounds<-c(
      Vv_wane_rate = 1/(365*1), 
      Vn_wane_rate = 1/(365*1), 
      rec_rate = 1/5,
      Vv_fail_rate=0, 
      Vv_symptom_rate=0, 
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_Av = 0.0, 
      beta_mod_An = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      Vv_0=0.0,
      Vn_0=0.0,
      Av_0=0.0,
      An_0=0.0,
      TDaP_waning_boost=0.0)
    
    upper.param.bounds<-c(
      Vv_wane_rate = 1/(365*50), 
      Vn_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      Vv_fail_rate=1, 
      Vv_symptom_rate=1, 
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_Av = 1, 
      beta_mod_An = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      Vv_0=1.0,
      Vn_0=1.0,
      Av_0=1.0,
      An_0=1.0,
      TDaP_waning_boost=5.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="Vn.equal.Vwp")
  {
    lower.param.bounds<-c(
      Vnwp_wane_rate = 1/(365*1), 
      Vap_wane_rate = 1/(365*1),
      rec_rate = 1/5,
      Vnwp_fail_rate=0, 
      Vnwp_symptom_rate=0, 
      Vap_fail_rate=0,
      Vap_symptom_rate=0, 
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_Anwp = 0.0, 
      beta_mod_Aap = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      Vnwp_0=0.0,
      Vap_0=0.0,
      Anwp_0=0.0,
      Aap_0=0.0)
    
    upper.param.bounds<-c(
      Vnwp_wane_rate = 1/(365*50), 
      Vap_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      Vnwp_fail_rate=1, 
      Vnwp_symptom_rate=1, 
      Vap_fail_rate=1, 
      Vap_symptom_rate=1, 
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_Anwp = 1, 
      beta_mod_Aap = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      Vnwp_0=1.0,
      Vap_0=1.0,
      Anwp_0=1.0,
      Aap_0=1.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="Vn.equal.Vwp.booster")
  {
    lower.param.bounds<-c(
      Vwp1_wane_rate = 1/(365*1), 
      Vwp2_wane_rate = 1/(365*1), 
      Vap_wane_rate = 1/(365*1),
      Vn_wane_rate = 1/(365*1), 
      rec_rate = 1/5,
      Vnwp_fail_rate=0, 
      Vnwp_symptom_rate=0, 
      Vap_fail_rate=0,
      Vap_symptom_rate=0, 
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_Anwp = 0.0, 
      beta_mod_Aap = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      Vwp1_0=0.0,
      Vwp2_0=0.0,
      Vap_0=0.0,
      Vn_0=0.0,
      Anwp_0=0.0,
      Aap_0=0.0,
      TDaP_waning_boost=0.0)
    
    upper.param.bounds<-c(
      Vwp1_wane_rate = 1/(365*50), 
      Vwp2_wane_rate = 1/(365*50), 
      Vap_wane_rate = 1/(365*50), 
      Vn_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      Vnwp_fail_rate=1, 
      Vnwp_symptom_rate=1, 
      Vap_fail_rate=1, 
      Vap_symptom_rate=1, 
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_Anwp = 1, 
      beta_mod_Aap = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      Vwp1_0=1.0,
      Vwp2_0=1.0,
      Vap_0=1.0,
      Vn_0=1.0,
      Anwp_0=1.0,
      Aap_0=1.0,
      TDaP_waning_boost=5.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="Vn.equal.Vap")
  {
    lower.param.bounds<-c(
      Vwp_wane_rate = 1/(365*1), 
      Vnap_wane_rate = 1/(365*1),
      rec_rate = 1/5,
      Vwp_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vnap_fail_rate=0,
      Vnap_symptom_rate=0, 
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_Awp = 0.0, 
      beta_mod_Anap = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      Vwp_0=0.0,
      Vnap_0=0.0,
      Awp_0=0.0,
      Anap_0=0.0)
    
    upper.param.bounds<-c(
      Vwp_wane_rate = 1/(365*50), 
      Vnap_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      Vwp_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vnap_fail_rate=1, 
      Vnap_symptom_rate=1, 
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_Awp = 1, 
      beta_mod_Anap = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      Vwp_0=1.0,
      Vnap_0=1.0,
      Awp_0=1.0,
      Anap_0=1.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="Vn.equal.Vap.booster")
  {
    lower.param.bounds<-c(
      Vwp1_wane_rate = 1/(365*1), 
      Vwp2_wane_rate = 1/(365*1), 
      Vap_wane_rate = 1/(365*1),
      Vn_wane_rate = 1/(365*1), 
      rec_rate = 1/5,
      Vwp_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vnap_fail_rate=0,
      Vnap_symptom_rate=0, 
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_Awp = 0.0, 
      beta_mod_Anap = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      Vwp1_0=0.0,
      Vwp2_0=0.0,
      Vap_0=0.0,
      Vn_0=0.0,
      Awp_0=0.0,
      Anap_0=0.0,
      TDaP_waning_boost=0.0)
    
    upper.param.bounds<-c(
      Vwp1_wane_rate = 1/(365*50), 
      Vwp2_wane_rate = 1/(365*50), 
      Vap_wane_rate = 1/(365*50), 
      Vn_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      Vwp_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vnap_fail_rate=1, 
      Vnap_symptom_rate=1, 
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_Awp = 1, 
      beta_mod_Anap = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      Vwp1_0=1.0,
      Vwp2_0=1.0,
      Vap_0=1.0,
      Vn_0=1.0,
      Awp_0=1.0,
      Anap_0=1.0,
      TDaP_waning_boost=5.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="none.equal")
  {
    lower.param.bounds<-c(
      Vwp_wane_rate = 1/(365*1), 
      Vap_wane_rate = 1/(365*1),
      Vn_wane_rate = 1/(365*1), 
      rec_rate = 1/5,
      Vwp_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vap_fail_rate=0,
      Vap_symptom_rate=0, 
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_Awp = 0.0, 
      beta_mod_Aap = 0.0, 
      beta_mod_An = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      Vwp_0=0.0,
      Vap_0=0.0,
      Vn_0=0.0,
      Awp_0=0.0,
      Aap_0=0.0,
      An_0=0.0)
    
    upper.param.bounds<-c(
      Vwp_wane_rate = 1/(365*50), 
      Vap_wane_rate = 1/(365*50), 
      Vn_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      Vwp_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vap_fail_rate=1, 
      Vap_symptom_rate=1, 
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_Awp = 1, 
      beta_mod_Aap = 1, 
      beta_mod_An = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      Vwp_0=1.0,
      Vap_0=1.0,
      Vn_0=1.0,
      Awp_0=1.0,
      Aap_0=1.0,
      An_0=1.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
  
  if (model=="none.equal.booster")
  {
    lower.param.bounds<-c(
      Vwp1_wane_rate = 1/(365*1), 
      Vwp2_wane_rate = 1/(365*1), 
      Vap_wane_rate = 1/(365*1),
      Vn_wane_rate = 1/(365*1), 
      rec_rate = 1/5,
      Vwp_fail_rate=0, 
      Vwp_symptom_rate=0, 
      Vap_fail_rate=0,
      Vap_symptom_rate=0, 
      Vn_fail_rate=0, 
      Vn_symptom_rate=0,
      beta0 = 1e-3, 
      beta1 = 0,
      beta_mod_Awp = 0.0, 
      beta_mod_Aap = 0.0, 
      beta_mod_An = 0.0, 
      rho=1e-9,
      sigmaSE=0,
      lag=0,
      S_0=0.0,
      I_0=0.0,
      Vwp1_0=0.0,
      Vwp2_0=0.0,
      Vap_0=0.0,
      Vn_0=0.0,
      Awp_0=0.0,
      Aap_0=0.0,
      An_0=0.0,
      TDaP_waning_boost=0.0)
    
    upper.param.bounds<-c(
      Vwp1_wane_rate = 1/(365*50), 
      Vwp2_wane_rate = 1/(365*50), 
      Vap_wane_rate = 1/(365*50), 
      Vn_wane_rate = 1/(365*50), 
      rec_rate = 1/31,
      Vwp_fail_rate=1, 
      Vwp_symptom_rate=1, 
      Vap_fail_rate=1, 
      Vap_symptom_rate=1, 
      Vn_fail_rate=1, 
      Vn_symptom_rate=1,
      beta0 = 1e4, 
      beta1 = 1,
      beta_mod_Awp = 1, 
      beta_mod_Aap = 1, 
      beta_mod_An = 1, 
      rho=1.0,
      sigmaSE=10,
      lag=1,
      S_0=1.0,
      I_0=1.0,
      Vwp1_0=1.0,
      Vwp2_0=1.0,
      Vap_0=1.0,
      Vn_0=1.0,
      Awp_0=1.0,
      Aap_0=1.0,
      An_0=1.0,
      TDaP_waning_boost=5.0)
    
    LHS<-sobolDesign(lower=lower.param.bounds,upper=upper.param.bounds,nseq=n.points)
    saveRDS(LHS,file=LHS.file.name)
  }
}
setwd(base.dir)

  
  