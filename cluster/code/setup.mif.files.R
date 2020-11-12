### setup mif runs  on cluster ###
library(pomp)

loc<-"US"
model<-"test.stoch"
subset.data<-"wP"
smooth.interval<-"two.weeks"

n.initial<-2000 #number of LHS samples
n.final<-250
jobs.per.node<-25 #number of LHS samples to analyze in the same script

### set directories for current analysis ###

job.name<-paste0(loc,".",model,".",subset.data,".",smooth.interval)
base.dir<-paste0("~/pertussis/",job.name)

### set directories for reading files to be re-written and run###
data.dir<-"~/pertussis/data" #directory containing data files
out.dir<-"~/pertussis/output" #directory for writing output files
lhs.dir<-"~/pertussis/lhs" #directory containing LHS files
code.dir<-"~/pertussis/code" #directory containing code files
batch.dir<-paste0(base.dir,"/batch")

### generate LHS if one doesn't already exist
setwd(code.dir)
source("lhs.gen.R")

setwd(base.dir)

if(!file.exists(paste0(job.name,".initial.q")))
{
  writeLines(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH -o init.pertussis.out",
      "#SBATCH -e init.pertussis.err",
      paste0("#SBATCH -J ",job.name,".init"), 
      paste0("#SBATCH --array=1-",n.initial/jobs.per.node),
      "#SBATCH --sockets-per-node=1",
      "#SBATCH --cores-per-socket=5",
      "#SBATCH -t 0-48:00:00",
      "#SBATCH --mail-type=END",
      "#SBATCH --mail-user=ifmiller@princeton.edu",
      "",
      paste0("cd ",batch.dir,"/dir.i$SLURM_ARRAY_TASK_ID"),
      "srun R CMD BATCH parameter.sweep.initial.R"
    ),
    con=paste0(job.name,".initial.q")
  )
}

setwd(base.dir)

if(!file.exists(paste0(job.name,".compile.initial.q")))
{
  writeLines(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH -o comp.init.pertussis.out",
      "#SBATCH -e comp.init.pertussis.err",
      paste0("#SBATCH -J ",job.name,".comp.init"), 
      "#SBATCH --array=1",
      "#SBATCH --sockets-per-node=1",
      "#SBATCH --cores-per-socket=1",
      "#SBATCH -t 0-00:05:00",
      "#SBATCH --mail-type=END",
      "#SBATCH --mail-user=ifmiller@princeton.edu",
      "",
      paste0("cd ",base.dir),
      "srun R CMD BATCH compile.initial.output.R"
    ),
    con=paste0(job.name,".comp.initial.q")
  )
}

if(!file.exists(paste0(job.name,".final.q")))
{
  writeLines(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH -o fin.pertussis.out",
      "#SBATCH -e fin.pertussis.err",
      paste0("#SBATCH -J ",job.name,".final"), 
      paste0("#SBATCH --array=1-",n.final/jobs.per.node),
      "#SBATCH --sockets-per-node=1",
      "#SBATCH --cores-per-socket=5",
      "#SBATCH -t 0-72:00:00",
      "#SBATCH --mail-type=END",
      "#SBATCH --mail-user=ifmiller@princeton.edu",
      "",
      paste0("cd ",batch.dir,"/dir.f$SLURM_ARRAY_TASK_ID"),
      "srun R CMD BATCH parameter.sweep.final.R"
    ),
    con=paste0(job.name,".final.q")
  )
}

if(!file.exists(paste0(job.name,".compile.final.q")))
{
  writeLines(
    c(
      "#!/bin/bash",
      "#",
      "#SBATCH -o comp.final.pertussis.out",
      "#SBATCH -e comp.final.pertussis.err",
      paste0("#SBATCH -J ",job.name,".comp.fin"), 
      "#SBATCH --array=1",
      "#SBATCH --sockets-per-node=1",
      "#SBATCH --cores-per-socket=1",
      "#SBATCH -t 0-00:05:00",
      "#SBATCH --mail-type=END",
      "#SBATCH --mail-user=ifmiller@princeton.edu",
      "",
      paste0("cd ",base.dir),
      "srun R CMD BATCH compile.final.output.R"
    ),
    con=paste0(job.name,".comp.final.q")
  )
}


### setup folders for batch job ###

if(!dir.exists(batch.dir)) {dir.create(batch.dir)}

for (i in 0:((n.initial/jobs.per.node)-1))
{
  setwd(code.dir)
  initial.param.sweep.lines<-readLines("parameter.sweep.initial.R")
  initial.param.sweep.lines[grep("n.initial<-",initial.param.sweep.lines)]<-paste0("n.initial<-",n.initial)
  initial.param.sweep.lines[grep("n.final<-",initial.param.sweep.lines)]<-paste0("n.final<-",n.final)
  initial.param.sweep.lines[grep("jobs.per.node<-",initial.param.sweep.lines)]<-paste0("jobs.per.node<-",jobs.per.node)
  initial.param.sweep.lines[grep("data.dir<-",initial.param.sweep.lines)]<-paste0("data.dir<-",'"',data.dir,'"')
  initial.param.sweep.lines[min(grep("out.dir<-",initial.param.sweep.lines))]<-paste0("out.dir<-",'"',out.dir,'"')
  initial.param.sweep.lines[grep("lhs.dir<-",initial.param.sweep.lines)]<-paste0("lhs.dir<-",'"',lhs.dir,'"')
  initial.param.sweep.lines[grep("code.dir<-",initial.param.sweep.lines)]<-paste0("code.dir<-",'"',code.dir,'"')
  initial.param.sweep.lines[grep("loc<-",initial.param.sweep.lines)]<-paste0("loc<-",'"',loc,'"')
  initial.param.sweep.lines[grep("model<-",initial.param.sweep.lines)]<-paste0("model<-",'"',model,'"')
  initial.param.sweep.lines[grep("subset.data<-",initial.param.sweep.lines)]<-paste0("subset.data<-",'"',subset.data,'"')
  initial.param.sweep.lines[grep("smooth.interval<-",initial.param.sweep.lines)]<-paste0("smooth.interval<-",'"',smooth.interval,'"')
  initial.param.sweep.lines[grep("start.job.index<-",initial.param.sweep.lines)]<-paste0("start.job.index<-",i*jobs.per.node+1)
  
  setwd(batch.dir)
  if(!dir.exists(paste0("dir.i",i+1))) {dir.create(paste0("dir.i",i+1))}
  setwd(paste0(getwd(),"/dir.i",i+1))
  writeLines(initial.param.sweep.lines,con="parameter.sweep.initial.R")
}

for (i in 0:((n.final/jobs.per.node)-1))
{
  setwd(code.dir)
  final.param.sweep.lines<-readLines("parameter.sweep.final.R")
  final.param.sweep.lines[grep("n.initial<-",final.param.sweep.lines)]<-paste0("n.initial<-",n.initial)
  final.param.sweep.lines[grep("n.final<-",final.param.sweep.lines)]<-paste0("n.final<-",n.final)
  final.param.sweep.lines[grep("jobs.per.node<-",final.param.sweep.lines)]<-paste0("jobs.per.node<-",jobs.per.node)
  final.param.sweep.lines[grep("data.dir<-",final.param.sweep.lines)]<-paste0("data.dir<-",'"',data.dir,'"')
  final.param.sweep.lines[min(grep("out.dir<-",final.param.sweep.lines))]<-paste0("out.dir<-",'"',out.dir,'"')
  final.param.sweep.lines[grep("lhs.dir<-",final.param.sweep.lines)]<-paste0("lhs.dir<-",'"',lhs.dir,'"')
  final.param.sweep.lines[grep("code.dir<-",final.param.sweep.lines)]<-paste0("code.dir<-",'"',code.dir,'"')
  final.param.sweep.lines[grep("loc<-",final.param.sweep.lines)]<-paste0("loc<-",'"',loc,'"')
  final.param.sweep.lines[grep("model<-",final.param.sweep.lines)]<-paste0("model<-",'"',model,'"')
  final.param.sweep.lines[grep("subset.data<-",final.param.sweep.lines)]<-paste0("subset.data<-",'"',subset.data,'"')
  final.param.sweep.lines[grep("smooth.interval<-",final.param.sweep.lines)]<-paste0("smooth.interval<-",'"',smooth.interval,'"')
  final.param.sweep.lines[grep("start.job.index<-",final.param.sweep.lines)]<-paste0("start.job.index<-",i*jobs.per.node+1)
  
  setwd(batch.dir)
  if(!dir.exists(paste0("dir.f",i+1))) {dir.create(paste0("dir.f",i+1))}
  setwd(paste0(getwd(),"/dir.f",i+1))
  writeLines(final.param.sweep.lines,con="parameter.sweep.final.R")
}

if(!dir.exists("~/pertussis/output/partial")) {dir.create("~/pertussis/output/partial")}
  
setwd(code.dir)
compile.initial.output.lines<-readLines("compile.initial.output.R")

compile.initial.output.lines[grep("n.initial<-",compile.initial.output.lines)]<-paste0("n.initial<-",n.initial)
compile.initial.output.lines[grep("n.final<-",compile.initial.output.lines)]<-paste0("n.final<-",n.final)
compile.initial.output.lines[grep("jobs.per.node<-",compile.initial.output.lines)]<-paste0("jobs.per.node<-",jobs.per.node)
compile.initial.output.lines[grep("data.dir<-",compile.initial.output.lines)]<-paste0("data.dir<-",'"',data.dir,'"')
compile.initial.output.lines[min(grep("out.dir<-",compile.initial.output.lines))]<-paste0("out.dir<-",'"',out.dir,'"')
compile.initial.output.lines[grep("lhs.dir<-",compile.initial.output.lines)]<-paste0("lhs.dir<-",'"',lhs.dir,'"')
compile.initial.output.lines[grep("code.dir<-",compile.initial.output.lines)]<-paste0("code.dir<-",'"',code.dir,'"')
compile.initial.output.lines[grep("loc<-",compile.initial.output.lines)]<-paste0("loc<-",'"',loc,'"')
compile.initial.output.lines[grep("model<-",compile.initial.output.lines)]<-paste0("model<-",'"',model,'"')
compile.initial.output.lines[grep("subset.data<-",compile.initial.output.lines)]<-paste0("subset.data<-",'"',subset.data,'"')
compile.initial.output.lines[grep("smooth.interval<-",compile.initial.output.lines)]<-paste0("smooth.interval<-",'"',smooth.interval,'"')
compile.initial.output.lines[grep("start.job.index<-",compile.initial.output.lines)]<-paste0("start.job.index<-",i*jobs.per.node+1)

setwd(base.dir)
writeLines(compile.initial.output.lines,con="compile.initial.output.R")

setwd(code.dir)
compile.final.output.lines<-readLines("compile.final.output.R")

compile.final.output.lines[grep("n.final<-",compile.final.output.lines)]<-paste0("n.final<-",n.final)
compile.final.output.lines[grep("n.final<-",compile.final.output.lines)]<-paste0("n.final<-",n.final)
compile.final.output.lines[grep("jobs.per.node<-",compile.final.output.lines)]<-paste0("jobs.per.node<-",jobs.per.node)
compile.final.output.lines[grep("data.dir<-",compile.final.output.lines)]<-paste0("data.dir<-",'"',data.dir,'"')
compile.final.output.lines[min(grep("out.dir<-",compile.final.output.lines))]<-paste0("out.dir<-",'"',out.dir,'"')
compile.final.output.lines[grep("lhs.dir<-",compile.final.output.lines)]<-paste0("lhs.dir<-",'"',lhs.dir,'"')
compile.final.output.lines[grep("code.dir<-",compile.final.output.lines)]<-paste0("code.dir<-",'"',code.dir,'"')
compile.final.output.lines[grep("loc<-",compile.final.output.lines)]<-paste0("loc<-",'"',loc,'"')
compile.final.output.lines[grep("model<-",compile.final.output.lines)]<-paste0("model<-",'"',model,'"')
compile.final.output.lines[grep("subset.data<-",compile.final.output.lines)]<-paste0("subset.data<-",'"',subset.data,'"')
compile.final.output.lines[grep("smooth.interval<-",compile.final.output.lines)]<-paste0("smooth.interval<-",'"',smooth.interval,'"')
compile.final.output.lines[grep("start.job.index<-",compile.final.output.lines)]<-paste0("start.job.index<-",i*jobs.per.node+1)

setwd(base.dir)
writeLines(compile.final.output.lines,con="compile.final.output.R")

