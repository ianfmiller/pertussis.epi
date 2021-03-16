loc<-"US" 
subset.data<-"all" 
smooth.interval<-"two.weeks" 
models<-c("waning.diff.protection.diff.symptoms.diff.transmission.diff",
          "waning.same.protection.diff.symptoms.diff.transmission.diff",
          "waning.diff.protection.same.symptoms.diff.transmission.diff",
          "waning.diff.protection.diff.symptoms.same.transmission.diff",
          "waning.diff.protection.diff.symptoms.diff.transmission.same",
          "waning.same.protection.same.symptoms.diff.transmission.diff",
          "waning.same.protection.diff.symptoms.same.transmission.diff",
          "waning.same.protection.diff.symptoms.diff.transmission.same",
          "waning.diff.protection.same.symptoms.same.transmission.diff",
          "waning.diff.protection.same.symptoms.diff.transmission.same",
          "waning.diff.protection.diff.symptoms.same.transmission.same",
          "waning.same.protection.same.symptoms.same.transmission.diff",
          "waning.same.protection.same.symptoms.diff.transmission.same",
          "waning.same.protection.diff.symptoms.same.transmission.same",
          "waning.diff.protection.same.symptoms.same.transmission.same",
          "waning.same.protection.same.symptoms.same.transmission.same")

n.params<-c(29,28,28,28,28,27,27,27,27,27,27,26,26,26,26,25)

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
