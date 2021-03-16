# define model
locs<-c("US")
models<-c(
  "waning.diff.protection.diff.symptoms.diff.transmission.diff",
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

model.files<-c("build.pomp.dddd.R",
               "build.pomp.sddd.R",
               "build.pomp.dsdd.R",
               "build.pomp.ddsd.R",
               "build.pomp.ddds.R",
               "build.pomp.ssdd.R",
               "build.pomp.sdsd.R",
               "build.pomp.sdds.R",
               "build.pomp.dssd.R",
               "build.pomp.dsds.R",
               "build.pomp.ddss.R",
               "build.pomp.sssd.R",
               "build.pomp.ssds.R",
               "build.pomp.sdss.R",
               "build.pomp.dsss.R",
               "build.pomp.ssss.R")
loc<-locs[1]
model<-models[1]
subset.data<-"all" 
smooth.interval<-"two.weeks" 
job.name<-paste0(loc,".",model,".",subset.data,".",smooth.interval)

# plot LHS
setwd("~/Documents/GitHub/pertussis.epi/results/LHS")
data<-readRDS(paste0(job.name,".LHS.RDS"))
initial.condition.cols<-grep("_0",colnames(data))
param.cols<-seq(1:dim(data)[2])[-c(1:7,initial.condition.cols)]
sub.param.cols<-param.cols[grep("beta",colnames(data)[param.cols])]
pairs(data[,sub.param.cols],col="grey",pch=16,cex=.4)
pairs(data[,initial.condition.cols],col="grey",pch=16,cex=.4)

# plot initial results
setwd("~/Documents/GitHub/pertussis.epi/results/initial.results")
data<-read.csv(paste0("initial.sweep.",job.name,".csv"))
initial.condition.cols<-grep("_0",colnames(data))
param.cols<-seq(1:dim(data)[2])[-c(1:7,initial.condition.cols)]
sub.param.cols<-param.cols[grep("beta",colnames(data)[param.cols])]
pairs(data[,sub.param.cols],col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],rgb(1,0,0,alpha=.5),rgb(.6,.6,.6,alpha=.5)),pch=16,cex=.4)
pairs(data[,initial.condition.cols],col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],rgb(1,0,0,alpha=.5),rgb(.6,.6,.6,alpha=.5)),pch=16,cex=.4)

#  plot mid results
setwd("~/Documents/GitHub/pertussis.epi/results/mid.results")
data<-read.csv(paste0("mid.sweep.",job.name,".csv"))
initial.condition.cols<-grep("_0",colnames(data))
param.cols<-seq(1:dim(data)[2])[-c(1:7,initial.condition.cols)]
sub.param.cols<-param.cols[grep("beta",colnames(data)[param.cols])]
pairs(data[,sub.param.cols],col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],rgb(1,0,0,alpha=.5),rgb(.6,.6,.6,alpha=.5)),pch=16,cex=.4)
pairs(data[,initial.condition.cols],col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],rgb(1,0,0,alpha=.5),rgb(.6,.6,.6,alpha=.5)),pch=16,cex=.4)

# plot final results
setwd("~/Documents/GitHub/pertussis.epi/results/final.results")
data<-read.csv(paste0("final.sweep.",job.name,".csv"))
initial.condition.cols<-grep("_0",colnames(data))
param.cols<-seq(1:dim(data)[2])[-c(1:7,initial.condition.cols)]
sub.param.cols<-param.cols[grep("beta",colnames(data)[param.cols])]
pairs(data[,sub.param.cols],col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],rgb(1,0,0,alpha=.5),rgb(.6,.6,.6,alpha=.5)),pch=16,cex=.4)
pairs(data[,initial.condition.cols],col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],rgb(1,0,0,alpha=.5),rgb(.6,.6,.6,alpha=.5)),pch=16,cex=.4)

# two parameter plot
## grey is LHS, red is initial run, blue is final run. Size proportional to likelihood
par(mfrow=c(1,1))
setwd("~/Documents/GitHub/pertussis.epi/results/LHS")
lhs.data<-readRDS(paste0(job.name,".LHS.RDS"))
plot(lhs.data$beta0,lhs.data$beta1,col="grey",pch=16)

setwd("~/Documents/GitHub/pertussis.epi/results/initial.results")
initial.data<-read.csv(paste0("initial.sweep.",job.name,".csv"))
initial.data<-subset(initial.data,loglik>=quantile(initial.data$loglik,.95))
points(initial.data$beta0,initial.data$beta1,col="blue",pch=16,cex=(initial.data$loglik-min(initial.data$loglik))/(max(initial.data$loglik)-min(initial.data$loglik)))

setwd("~/Documents/GitHub/pertussis.epi/results/mid.results")
mid.data<-read.csv(paste0("mid.sweep.",job.name,".csv"))
mid.data<-subset(mid.data,loglik>=quantile(mid.data$loglik,.95))
points(mid.data$beta0,mid.data$beta1,col="purple",pch=16,cex=(mid.data$loglik-min(mid.data$loglik))/(max(mid.data$loglik)-min(mid.data$loglik)))

setwd("~/Documents/GitHub/pertussis.epi/results/final.results")
final.data<-read.csv(paste0("final.sweep.",job.name,".csv"))
final.data<-subset(final.data,loglik>=quantile(final.data$loglik,.95))
points(final.data$beta0,final.data$beta1,col="red",pch=16,cex=1)

# plot simulations

## setup
setwd("~/Documents/GitHub/pertussis.epi/results/LHS")
LHS<-readRDS(paste0(job.name,".LHS.RDS"))
data.dir<-"~/Documents/GitHub/pertussis.epi/cluster/data"
setwd("~/Documents/GitHub/pertussis.epi/cluster/code")
source("prep.data.covar.R")
setwd("~/Documents/GitHub/pertussis.epi/cluster/code/models")
source(model.files[which(models==model)])

setwd("~/Documents/GitHub/pertussis.epi/results/final.results")
fit.params<-read.csv(paste0("final.sweep.",job.name,".csv"))
initial.condition.cols<-grep("_0",colnames(fit.params))
param.cols<-seq(1:dim(fit.params)[2])[-c(1:7,initial.condition.cols)]
params<-fit.params[order(-fit.params$loglik)[4],c(param.cols,initial.condition.cols)]

## plot individuals simulations

simulate(m1,params=params,nsim=8,format="d",include.data=TRUE,verbose=T)->out


par(mfrow=c(3,3),mar=c(3,3,1,0.5),mgp=c(1.5,.5,0))
plot(subset(out,.id=="data")$time,subset(out,.id=="data")$incidence,type="l",xlab="time",ylab="incidence",main="data",cex.axis=.5)
for (i in 1:8)
{
  plot(subset(out,.id==i)$time,subset(out,.id==i)$incidence,type="l",xlab="time",ylab="incidence",main=paste("sim = ",i,sep=""),cex.axis=.5,col="blue")
}

## plot bounds
par(mfrow=c(1,1))
simulate(m1,params=params,nsim=1000,format="d",include.data=TRUE)->out
out<-out[,c("time","incidence")]
upper<-c()
lower<-c()
for(t in dates)
{
  out.sub<-subset(out,time==t)
  bounds<-quantile(out.sub$incidence,probs=c(.05,.95))
  lower<-c(lower,bounds[1])
  upper<-c(upper,bounds[2])
}

plot(m1@times,c(m1@data),type="l",xlab="time",ylab="incidence",main="data",cex.axis=.5,ylim=c(0,max(upper)))
points(m1@times,upper,type="l",col="red")
points(m1@times,lower,type="l",col="red")
polygon(c(m1@times,rev(m1@times)),c(lower,rev(upper)),col="pink",border = NA)
points(m1@times,c(m1@data),type="l")

## visualize individual mif
i<-395 ### index of mif of interest
setwd("~/Downloads/results.out.2") ### directory containing mif files
file.name<-paste0(model,".",loc,".",subset.data,".",smooth.interval,".iter",i,".final.mif.RDS")
mod<-readRDS(file.name)
plot(1:length(mod$mif@traces[,"loglik"]),mod$mif@traces[,"beta_mod_An"],type="l")

