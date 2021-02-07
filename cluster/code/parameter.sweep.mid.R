#### final LHS parameter sweep--stochastic model ###
library(doRNG)
library(doParallel)
library(pomp)
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
smooth.interval<-"four.weeks" #time period over which to smooth data, options are "two.weeks" "four.weeks" and "none"
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
if(model=="none.equal.dev") {source("build.pomp.waning.diff.protection.diff.symptoms.diff.transmission.diff.R")}

### load start points from initial sweep

initial.sweep<-read.csv(paste0(out.dir,"/initial.sweep.",loc,".",model,".",subset.data,".",smooth.interval,".csv"))
initial.sweep<-initial.sweep[order(initial.sweep$loglik,decreasing = T),]
initial.sweep<-initial.sweep[1:n.mid,]
job.set<-initial.sweep$lhs.row

ncores=detectCores()
registerDoParallel(cores=ncores)

foreach(i=0:(jobs.per.node-1), .inorder=F, .combine = "rbind") %dorng% {
  print(paste("starting i =",i+start.job.index))
  
  job.index<-start.job.index+i
  lhs.samp<-job.set[job.index]
  setwd(out.dir)
  m1<-readRDS(paste0(model,".",loc,".",subset.data,".",smooth.interval,".iter",lhs.samp,".initial.mif.RDS"))
  m1<-m1$mif
  continue(m1,Nmif=15)->m2
  print(paste0("i = ",i+start.job.index," mif complete; time = ",Sys.time()))
  ll <- replicate(n=10,logLik(pfilter(m2,Np=10000)))
  print(paste("finished i =",i+start.job.index,"; time = ",Sys.time()))
  m2<-list(mif=m2,ll=logmeanexp(ll,se=TRUE))
  setwd(out.dir)
  saveRDS(m2,file=paste(model,loc,subset.data,smooth.interval,paste("iter",lhs.samp,sep=""),"mid.mif.RDS",sep="."))
  data.frame("loc"=loc,"model"=model,"subset.data"=subset.data,"smooth.interval"=smooth.interval,"lhs.row"=lhs.samp,"loglik"=m2$ll[[1]],"se"=m2$ll[[2]],rbind(coef(m2$mif)))
}->mf.out

setwd(out.dir)
if(!dir.exists("partial")) {dir.create("partial")}
setwd(partial.out.dir)
out.name<-paste0(start.job.index,"mid.sweep.",loc,".",model,".",subset.data,".",smooth.interval,".csv")
if (file.exists(out.name)) {warning("old file overwritten")}
write.csv(mf.out,file=out.name,row.names = F)
