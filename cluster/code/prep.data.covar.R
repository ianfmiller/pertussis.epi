#### build pertussis incidence and covariates datasets ###
library(lubridate)
setwd(data.dir)

### load data files ###

data<-read.csv("raw.incidence.csv",row.names = 1)
yearly.totals<-read.csv("yearly totals.csv",row.names = 1,check.names=FALSE)
demog<-read.csv("demog.csv")
vacc<-read.csv("vacc.DTP3.rate.at.19.to.35.months.csv",row.names = 1,check.names = F)

if(loc=="multi.test") {orig.loc<-loc; loc<-c("Virginia","Maryland","District.of.Columbia")}

### formal data ###

data<-t(data)
data<-gsub(",","",data)
data<-gsub(" ","",data)

data[which(data=="")]<-0
data[which(is.na(data))]<-0
data[which(data=="U")]<-0 
data[which(data<0)]<-0
dates<-as.Date(gsub("X","",rownames(data)),tryFormats=c("%m.%d.%Y"))

### correct reported incidence by yearly total ###
data<-data.frame(apply(data, 2, function(x) as.numeric(as.character(x))))
data<-as.data.frame(data[,loc])
colnames(data)<-loc
rownames(data)<-as.character(dates)

new.data<-c()
for (k in 1:dim(data)[2])
{
  sub.data<-as.numeric(data[,k])
  new.sub.data<-c()
  for (i in 1:length(sub.data))
  {
    sub.loc<-colnames(data)[k]
    reported.data<-as.numeric(as.character(sub.data[i]))
    sum.year.reported.data<-sum(as.numeric(as.character(sub.data[which(year(dates)==year(dates[i]))])))
    year.total<-yearly.totals[as.character(year(dates[i])),sub.loc]  
    new.sub.data<-c(new.sub.data,reported.data/sum.year.reported.data*year.total)
    if(sum.year.reported.data==0) {new.sub.data[length(new.sub.data)]<-0}
  }
  new.data<-cbind(new.data,new.sub.data)
}

colnames(new.data)<-colnames(data)
rownames(new.data)<-rownames(data)
data<-new.data

if(dim(data)[2]>1)
{
  new.data<-as.data.frame(rowSums(data))
  colnames(new.data)<-orig.loc
  data<-new.data
}

### smooth data ###

smooth.data<-function(smooth.interval,x)
{
  if(smooth.interval=="two.weeks")
  { 
    dates<-as.Date(rownames(x),tryFormats=c("%Y-%m-%d"))
    if (dim(x)[1] %% 2 == 1) {x<-x[-c(dim(x)[1]),];dates<-dates[-length(dates)]}
    dates<-dates[seq(1,length(dates),2)]
    i<-1
    out.data<-c()
    while (i<=length(x))
    {
      out.data<-rbind(out.data,x[i]+x[i+1])
      i<-i+2
    }
  }
  
  if(smooth.interval=="four.weeks")
  { 
    dates<-as.Date(rownames(x),tryFormats=c("%Y-%m-%d"))
    if (!dim(x)[1] %% 4 == 0) {x<-x[-((dim(x)[1]):(dim(x)[1]-(dim(x)[1] %% 4)+1)),];dates<-dates[-(length(dates):(length(dates)-length(dates)%%4))] }
    dates<-dates[seq(1,length(dates),4)]
    i<-1
    out.data<-c()
    while (i<=length(x))
    {
      out.data<-rbind(out.data,x[i]+x[i+1]+x[i+2]+x[i+3])
      i<-i+4
    }
  }
  
  if(smooth.interval=="twenty.weeks") ## for tests
  { 
    dates<-as.Date(rownames(x),tryFormats=c("%Y-%m-%d"))
    if (!dim(x)[1] %% 20 == 0) {x<-x[-((dim(x)[1]):(dim(x)[1]-(dim(x)[1] %% 20)+1)),];dates<-dates[-(length(dates):(length(dates)-length(dates)%%20))] }
    dates<-dates[seq(1,length(dates),20)]
    i<-1
    out.data<-c()
    while (i<=length(x))
    {
      out.data<-rbind(out.data,x[i]+x[i+1]+x[i+2]+x[i+3]+x[i+4]+x[i+5]+x[i+6]+x[i+7]+x[i+8]+x[i+9]+x[i+10]+x[i+11]+x[i+12]+x[i+13]+x[i+14]+x[i+15]+x[i+16]+x[i+17]+x[i+18]+x[i+19]+x[i+20])
      i<-i+20
    }
  }

  rownames(out.data)<-as.character(dates)
  colnames(out.data)<-ifelse(length(loc)>1,orig.loc,loc)
  return(list(out.data,dates))
}

smoothed.data.obj<-smooth.data(smooth.interval,data)
data<-smoothed.data.obj[[1]]
dates<-smoothed.data.obj[[2]]

### subset data by period ###

if(subset.data=="all")
  {  
    dates<-as.numeric(dates)
    data<-data[1:length(data)]
  }

if(subset.data=="wP")
  {  
    date.ceiling<-max(which(dates<as.Date("1997-01-01")))
    dates<-dates[1:date.ceiling]
    data<-data[1:date.ceiling]
    dates<-as.numeric(dates)
  }
  
if(subset.data=="aP")
  {  
    date.floor<-min(which(dates>=as.Date("1997-01-01")))
    dates<-dates[date.floor:length(dates)]
    data<-data[date.floor:length(dates)]
    dates<-as.numeric(dates)
  }


### make final data object ###

data<-data.frame(time=dates,incidence=round(data))

### prep demography data ###

demog<-subset(demog,State %in% loc)

if(length(loc)>1)
{
  new.Year<-unique(demog$Year)
  new.State<-rep(orig.loc,times=length(new.Year))
  new.Births<-aggregate(demog$Births,by=list(Year=demog$Year),FUN=sum)$x
  new.Deaths<-aggregate(demog$Deaths,by=list(Year=demog$Year),FUN=sum)$x
  new.Pop<-aggregate(demog$Pop,by=list(Year=demog$Year),FUN=sum)$x
  demog<-data.frame(Year=new.Year,State=new.State,Births=new.Births,Deaths=new.Deaths,Pop=new.Pop)
}

years<-as.Date(paste(demog$Year,"-01-01",sep="")) #convert year to date object
years<-as.numeric(years) #convert years to numeric

### prep covaraiates ###

birth.rate<-predict(smooth.spline(years,demog$Births/demog$Pop,df=length(years)),x=dates)$y

death.rate<-predict(smooth.spline(years,demog$Deaths/demog$Pop,df=length(years)),x=dates)$y

pop<-predict(smooth.spline(years,demog$Pop,df=length(years)),x=dates)$y

## prep vaccination rate data ##

# data for US is mostly complete since 1982
# Data for states/regions is incomplete from 1982-1994, but complete from 1995-2017. 
# Rates from 1982-1994 are calculated as the U.S. rate +/- the difference between the
# U.S. rate and the state/mean region in 1995.

  vacc.years<-as.Date(paste(colnames(vacc),"-01-01",sep="")) #convert year to date object
  vacc.years<-vacc.years[-which(is.na(vacc["US",]))] #cut years with missing data
  vacc.years<-as.numeric(vacc.years) #convert years to numeric
  US.rates<-as.numeric(vacc["US",])[-which(is.na(vacc["US",]))] #cut missing data
  vacc.data<-cbind(vacc.years,US.rates) #create vacc rate object
  colnames(vacc.data)<-c("year","US.rate")
  vacc.data<-as.data.frame(vacc.data)
  vacc.data$US.rate<-vacc.data$US.rate/100

if(all(loc=="US"))
{
  vacc.rate<-predict(smooth.spline(vacc.data$year,vacc.data$US.rate),x=dates)$y
}
  
if(!all(loc=="US"))
{
  loc.rates<-colSums(vacc[loc,])/length(loc)
  loc.rates<-loc.rates[-which(is.na(vacc["US",]))]
  vacc.data<-cbind(vacc.data,loc.rates/100)
  colnames(vacc.data)[3]<-"loc.rate"
  
  diff<-vacc.data[9,2]-vacc.data[9,3]
  vacc.data[1:8,"loc.rate"]<-vacc.data[1:8,"US.rate"]+diff
  
  vacc.rate<-predict(smooth.spline(vacc.data$year,vacc.data$loc.rate),x=dates)$y
}


### merge covariates

covartable<-data.frame(time=dates,birth_rate=birth.rate,death_rate=death.rate,pop=pop,vacc_rate=vacc.rate)

### restore original loc value if needed

if(exists("orig.loc")) {loc<-orig.loc}

