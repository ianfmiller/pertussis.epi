rm(list=ls())

setwd("~/pertussis/code")
analysis.setup.lines<-readLines("setup.mif.files.R")

locs<-c("US")
models<-c("all.equal","none.equal","Vn.equal.Vwp")
subset.data<-"all" #vaccine era to subset. options are "all" "wP" and "aP"
smooth.interval<-"two.weeks"  #time period over which to smooth data, options are "two.weeks" "four.weeks" and "none"
for(model in models)
{
for(loc in locs) #add nested loops if needed for other models, subsets, etc.
  {
    analysis.setup.lines[min(grep('loc<-',analysis.setup.lines))]<-paste0("loc<-",'"',loc,'"')
    analysis.setup.lines[min(grep('model<-',analysis.setup.lines))]<-paste0("model<-",'"',model,'"')
    analysis.setup.lines[min(grep('subset.data<-',analysis.setup.lines))]<-paste0("subset.data<-",'"',subset.data,'"')
    analysis.setup.lines[min(grep('smooth.interval<-',analysis.setup.lines))]<-paste0("smooth.interval<-",'"',smooth.interval,'"')
  
    job.name<-paste0(loc,".",model,".",subset.data,".",smooth.interval)
    setwd("~/pertussis")
    if (!dir.exists(job.name)) {dir.create(job.name)}
    setwd(job.name)
    writeLines(analysis.setup.lines,con="setup.mif.files.R")
    source("setup.mif.files.R")
    setwd(paste0("~/pertussis/",job.name))
    
    cmd<-paste0("jid1=$(sbatch ",paste0(job.name,".initial.q);"),
                " jid2=$(sbatch --dependency=afterany:${jid1##* } ",paste0(job.name,".comp.initial.q);"),
                " jid3=$(sbatch --dependency=afterany:${jid2##* } ",paste0(job.name,".mid.q);"),
                " jid4=$(sbatch --dependency=afterany:${jid3##* } ",paste0(job.name,".comp.mid.q);"),
                " jid5=$(sbatch --dependency=afterany:${jid4##* } ",paste0(job.name,".final.q);"),
                " sbatch --dependency=afterany:${jid5##* } ",paste0(job.name,".comp.final.q"))
    system(cmd)
}
}

