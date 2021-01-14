loc<-"US" 
subset.data<-"all" 
smooth.interval<-"two.weeks" 
models<-c("none.equal","all.equal","Vn.equal.Vwp")
n.params<-c(26,14,20)

logliks<-c()
for(model in models)
{
  job.name<-paste0(loc,".",model,".",subset.data,".",smooth.interval)
  fit.params<-read.csv(paste0("~/Documents/GitHub/pertussis.epi/results/final.results/final.sweep.",job.name,".csv"))
  loglik<-fit.params[order(-fit.params$loglik)[1],"loglik"]
  logliks<-c(logliks,loglik)
}

AICs<-2*n.params-2*logliks
names(AICs)<-models
AICs
