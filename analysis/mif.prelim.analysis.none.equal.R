loc<-"US" 
model<-"none.equal"
subset.data<-"all" 
smooth.interval<-"two.weeks" 
job.name<-paste0(loc,".",model,".",subset.data,".",smooth.interval)

## plot LHS
setwd("~/Dropbox (Princeton)/USPertussis/pomp/results/LHS")
data<-readRDS(paste0(job.name,".LHS.RDS"))
pairs(~beta0+beta1+beta_mod_Awp+rho+beta_mod_Aap+beta_mod_An+rec_rate,data=data,col="grey",pch=16,cex=.4)
pairs(~S_0+I_0+Vwp_0+Vap_0+Vn_0+Awp_0+Aap_0+An_0,data=data,col="grey",pch=16,cex=.4)

## plot initial results
setwd("~/Dropbox (Princeton)/USPertussis/pomp/results/initial results")
data<-read.csv(paste0("initial.sweep.",job.name,".csv"))
pairs(~beta0+beta1+beta_mod_Awp+rho+beta_mod_Aap+beta_mod_An+rec_rate,data=data,col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],rgb(1,0,0,alpha=.5),rgb(.6,.6,.6,alpha=.5)),pch=16,cex=.4)
pairs(~S_0+I_0+Vwp_0+Vap_0+Vn_0+Awp_0+Aap_0+An_0,data=data,col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],rgb(1,0,0,alpha=.5),rgb(.6,.6,.6,alpha=.5)),pch=16,cex=.4)
     
## plot final results
setwd("~/Dropbox (Princeton)/USPertussis/pomp/results/final results")
data<-read.csv(paste0("final.sweep.",job.name,".csv"))
pairs(~beta0+beta1+beta_mod_Awp+rho+beta_mod_Aap+beta_mod_An+rec_rate,data=data,col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],"red","grey"),pch=16,cex=.4)
pairs(~S_0+I_0+Vwp_0+Vap_0+Vn_0+Awp_0+Aap_0+An_0,data=data,col=ifelse((data$loglik) > quantile(data$loglik,.9)[[1]],"red","grey"),pch=16,cex=.4)

#### two parameter plot
## grey is LHS, red is initial run, blue is final run. Size proportional to likelihood
par(mfrow=c(1,1))
setwd("~/Dropbox (Princeton)/USPertussis/pomp/results/LHS")
data<-readRDS(paste0(job.name,".LHS.RDS"))
plot(data$beta0,data$beta1,col="grey",pch=16)

setwd("~/Dropbox (Princeton)/USPertussis/pomp/results/initial results")
data<-read.csv(paste0("initial.sweep.",job.name,".csv"))
data<-subset(data,loglik>=quantile(data$loglik,.95))
points(data$beta0,data$beta1,col="red",pch=16,cex=(data$loglik-min(data$loglik))/(max(data$loglik)-min(data$loglik)))


setwd("~/Dropbox (Princeton)/USPertussis/pomp/results/final results")
data<-read.csv(paste0("final.sweep.",job.name,".csv"))
data<-subset(data,loglik>=quantile(data$loglik,.95))
points(data$beta0,data$beta1,col="blue",pch=16,cex=(data$loglik-min(data$loglik))/(max(data$loglik)-min(data$loglik)))

### plot simulations
setwd("~/Dropbox (Princeton)/USPertussis/pomp/results/LHS")
LHS<-readRDS(paste0(job.name,".LHS.RDS"))
data.dir<-"~/Dropbox (Princeton)/USPertussis/pomp/cluster/data"
setwd("~/Dropbox (Princeton)/USPertussis/pomp/cluster/code")
source("prep.data.covar.R")
setwd("~/Downloads")
source("~/Dropbox (Princeton)/USPertussis/pomp/cluster/code/build.pomp.none.equal.R")

setwd("~/Dropbox (Princeton)/USPertussis/pomp/results/final results")
data<-read.csv(paste0("final.sweep.",job.name,".csv"))
params<-data[order(-data$loglik)[1],8:33]

simulate(m1,params=params,nsim=8,format="d",include.data=TRUE,verbose=T)->out


par(mfrow=c(3,3),mar=c(3,3,1,0.5),mgp=c(1.5,.5,0))
plot(subset(out,.id=="data")$time,subset(out,.id=="data")$incidence,type="l",xlab="time",ylab="incidence",main="data",cex.axis=.5)
for (i in 1:8)
{
  plot(subset(out,.id==i)$time,subset(out,.id==i)$incidence,type="l",xlab="time",ylab="incidence",main=paste("sim = ",i,sep=""),cex.axis=.5,col="blue")
}

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

i<-1565
setwd("~/Downloads/tout")
file.name<-paste0(model,".",loc,".",subset.data,".",smooth.interval,".iter",i,"..final.mif.RDS")
mod<-readRDS(file.name)
plot(0:950,mod$mif@traces[,"beta0"])
abline(v=350)

