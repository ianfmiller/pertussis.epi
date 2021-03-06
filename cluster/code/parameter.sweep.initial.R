#### LHS parameter sweep--stochastic model ###
library(doRNG)
library(doParallel)
library(R.utils)
set.seed(13548996)
### set job iteration ###

n.initial<-2000 #number of LHS samples
n.mid<-200  #number of LHS samples
n.final<-20 #number of LHS samples
jobs.per.node<-20 #number of LHS samples to analyze in the same script
start.job.index<-1 #first index parameter set to analyze

### set region, smoothing window, vaccine era, model
loc<-"US" #loc is the region or vector of regions to analyze NEED TO MAKE SURE THIS WORKS FOR DC, NY, AND NYC
model<-"test.stoch"
subset.data<-"wP" #vaccine era to subset. options are "all" "wP" and "aP"
smooth.interval<-"two.weeks" #time period over which to smooth data, options are "two.weeks" "four.weeks" and "none"
job.name<-paste0(loc,".",model,".",subset.data,".",smooth.interval)

### set directories for reading files to be re-written and run###
data.dir<-"~/pertussis/data" #directory containing data files
out.dir<-"~/pertussis/output" #directory for writing output files
partial.out.dir<-"~/pertussis/output/partial" #directory for writing partial output files
lhs.dir<-"~/pertussis/lhs" #directory containing LHS files
code.dir<-"~/pertussis/code" #directory containing code files
base.dir<-paste0("~/pertussis/",job.name)
batch.dir<-paste0(base.dir,"/batch")

### create LHS parameter sets
setwd(code.dir)
source("lhs.gen.R")

### prep data, load model
setwd(code.dir)
source("prep.data.covar.R")
setwd(code.dir)
setwd("models")
if(model=="waning.diff.protection.diff.symptoms.diff.transmission.diff") {source("build.pomp.dddd.R")}

if(model=="waning.same.protection.diff.symptoms.diff.transmission.diff") {source("build.pomp.sddd.R")}
if(model=="waning.diff.protection.same.symptoms.diff.transmission.diff") {source("build.pomp.dsdd.R")}
if(model=="waning.diff.protection.diff.symptoms.same.transmission.diff") {source("build.pomp.ddsd.R")}
if(model=="waning.diff.protection.diff.symptoms.diff.transmission.same") {source("build.pomp.ddds.R")}

if(model=="waning.same.protection.same.symptoms.diff.transmission.diff") {source("build.pomp.ssdd.R")}
if(model=="waning.same.protection.diff.symptoms.same.transmission.diff") {source("build.pomp.sdsd.R")}
if(model=="waning.same.protection.diff.symptoms.diff.transmission.same") {source("build.pomp.sdds.R")}
if(model=="waning.diff.protection.same.symptoms.same.transmission.diff") {source("build.pomp.dssd.R")}
if(model=="waning.diff.protection.same.symptoms.diff.transmission.same") {source("build.pomp.dsds.R")}
if(model=="waning.diff.protection.diff.symptoms.same.transmission.same") {source("build.pomp.ddss.R")}

if(model=="waning.same.protection.same.symptoms.same.transmission.diff") {source("build.pomp.sssd.R")}
if(model=="waning.same.protection.same.symptoms.diff.transmission.same") {source("build.pomp.ssds.R")}
if(model=="waning.same.protection.diff.symptoms.same.transmission.same") {source("build.pomp.sdss.R")}
if(model=="waning.diff.protection.same.symptoms.same.transmission.same") {source("build.pomp.dsss.R")}

if(model=="waning.same.protection.same.symptoms.same.transmission.same") {source("build.pomp.ssss.R")}


### analyze LHS sample ###

ncores=detectCores()
registerDoParallel(cores=ncores)

foreach(i=0:(jobs.per.node-1), .inorder=F, .combine = "rbind") %dorng% {
  print(paste("starting i =",i+start.job.index,"; time = ",Sys.time()))
  job.index<-start.job.index+i
  params<-params.mat[job.index,]
  #m2<-mif2(m1,Nmif=25,params=params,rw.sd=rw.sd,cooling.fraction.50=0.25, Np=1000,cooling.type="geometric") # for two.week smooth interval
  m2<-mif2(m1,Nmif=10,params=params,rw.sd=rw.sd,cooling.fraction.50=0.5, Np=1000,cooling.type="geometric") # for four.week smooth interval
  print(paste0("i = ",i+start.job.index," mif complete; time = ",Sys.time()))
  ### get liklihood with check to kill slow calculation of extremely low liklihood 
  ll <- tryCatch(withTimeout(replicate(n=10,logLik(pfilter(m2,Np=10000))),timeout = 60*60*20),error=function(e) {ll.tmp<- rep(-888e10,times=10)})
  print(paste("finished i =",i+start.job.index,"; time = ",Sys.time()))
  m2<-list(mif=m2,ll=logmeanexp(ll,se=TRUE))
  setwd(out.dir) 
  saveRDS(m2,file=paste(model,loc,subset.data,smooth.interval,paste("iter",job.index,sep=""),"initial.mif.RDS",sep="."))
  data.frame("loc"=loc,"model"=model,"subset.data"=subset.data,"smooth.interval"=smooth.interval,"lhs.row"=i+start.job.index,"loglik"=m2$ll[[1]],"se"=m2$ll[[2]],rbind(coef(m2$mif)))
}->mf.out

setwd(out.dir)
if(!dir.exists("partial")) {dir.create("partial")}
setwd(partial.out.dir)
out.name<-paste0(start.job.index,"initial.sweep.",loc,".",model,".",subset.data,".",smooth.interval,".csv")
if (file.exists(out.name)) {warning("old file overwritten")}
write.csv(mf.out,file=out.name,row.names = F)
