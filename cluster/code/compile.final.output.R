### compile data from final run ###
## data written in separate files for each core to avoid errors associated with multiple cores writing to the same file simultaneously

n.initial<-2000 #number of LHS samples
n.final<-250 #n best samples to analyze further
jobs.per.node<-25 #number of LHS samples to analyze in the same script


### set region, smoothing window, vaccine era, model
loc<-"US" #loc is the region or vector of regions to analyze NEED TO MAKE SURE THIS WORKS FOR DC, NY, AND NYC
model<-"test.stoch"
subset.data<-"wP" #vaccine era to subset. options are "all" "wP" and "aP"
smooth.interval<-"two.weeks" #time period over which to smooth data, options are "two.weeks" "four.weeks" and "none"
job.name<-paste0(loc,".",model,".",subset.data,".",smooth.interval)

### set directories for reading files to be re-written and run###
data.dir<-"~/pertussis/data" #directory containing data files
out.dir<-"~/pertussis/output" #directory for writing output files
partial.out.dir<-paste0(out.dir,"/partial") #directory for writing partial output files
lhs.dir<-"~/pertussis/lhs" #directory containing LHS files
code.dir<-"~/pertussis/code" #directory containing code files
base.dir<-paste0("~/pertussis/",job.name)
batch.dir<-paste0(base.dir,"/batch")


setwd(partial.out.dir)
output<-c()
for (i in seq(1,n.final,jobs.per.node))
{
  output.chunk<-read.csv(paste0(i,"final.sweep.",loc,".",model,".",subset.data,".",smooth.interval,".csv"))
  output<-rbind(output,output.chunk)
}

setwd(out.dir)
write.csv(output,file=paste0("final.sweep.",loc,".",model,".",subset.data,".",smooth.interval,".csv"),row.names = F)

setwd(partial.out.dir)
for (i in seq(1,n.final,jobs.per.node))
{
  file.remove(paste0(i,"final.sweep.",loc,".",model,".",subset.data,".",smooth.interval,".csv"))
}


