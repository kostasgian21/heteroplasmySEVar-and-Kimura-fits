# Current code, uncomment for parameters
# !/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 3) {
  stop("Needing 3 arguments")
}

dataType = myargs[1]
numberOfMCMCsamples=as.numeric(myargs[2])
nSamplesNoisedData=as.numeric(myargs[3])


library(kimura)
library(heteroplasmy)
library(foreach)




# #comment bellow for not-parallel code (and also line 151). Change also %dopar% to %do% in line 124
# library(doParallel)
# library(parallel)
# cores=4
# cl=makeCluster(cores)
# 
# registerDoParallel(cl)
# 
# clusterEvalQ(cl, expr={library(kimura)
#   library(heteroplasmy)})





# #initial parameters
# dataType="Offsprings" # possible values PGC, Oocytes, Offsprings
# numberOfMCMCsamples=50
# nSamplesNoisedData=2

mouseLine="freyer"

readFileName=paste(mouseLine,dataType,".txt",sep="")
mouseData <- read.table(readFileName, sep = "\t", header = T)

if (dataType=="Offsprings" || dataType=="PGC") {
  mouseDataHeaders=read.table(readFileName, sep = "\t", 
                              header = F,nrows=1)
  mouseDataHeaders1=which(mouseDataHeaders<=0.6)
  mouseDataHeaders2=which(mouseDataHeaders<=0.7 & mouseDataHeaders>0.6)
  mouseDataHeaders3=which(mouseDataHeaders>0.7)
  
  mouseData1=vector()
  for (i in mouseDataHeaders1) {
    mouseData1=c(mouseData1,mouseData[which(!is.na(mouseData[,i])),i])
  }
  
  mouseData2=vector()
  for (i in mouseDataHeaders2) {
    mouseData2=c(mouseData2,mouseData[which(!is.na(mouseData[,i])),i])
  }
  
  mouseData3=vector()
  for (i in mouseDataHeaders3) {
    mouseData3=c(mouseData3,mouseData[which(!is.na(mouseData[,i])),i])
  }
  numData=3
}else if(dataType=="Oocytes"){
  numData=nrow(mouseData)
}


pvals= data.frame(Moments=numeric(),
                  maxlik=numeric())

pvalsRange=data.frame(Min=numeric(),
                      Max=numeric())

for (z in 1:numData) {
  if (dataType=="Offsprings" || dataType=="PGC") {
  rhs=paste("X.1=","mouseData",z,sep="")
  eval(parse(text=rhs))
  }else{
    X.1=mouseData[which(!is.na(mouseData[,z])),z]
  }
  n=length(X.1)
  
  
  X.1[which(X.1>1)]=1
  X.1[which(X.1<0)]=0
  pvalstmp= data.frame(Moments=numeric(),
                       maxlik=numeric())
  mm <- foreach (i=1:nSamplesNoisedData, .combine=rbind) %do% {
    aa=kimura::test_kimura(X.1,num_MC = numberOfMCMCsamples)
    bb=test_kimura_par(X.1,p=estimate_parameters_ks(X.1)[1],b=estimate_parameters_ks(X.1)[2],num_MC = numberOfMCMCsamples)
    c(aa$p.value,bb$p.value)
  }
  pvalstmp[1:nSamplesNoisedData,]=mm[1:nSamplesNoisedData,]
  pvals[z,]=c(mean(pvalstmp[,1]),mean(pvalstmp[,2]))
  
  # range of MoM pvals
  pvalsRange[z,]=c(min(pvalstmp[,1]),max(pvalstmp[,1]))
  write.table(pvals, paste(dataType,"_",mouseLine,"pvals.txt",sep = ""), sep = " ",row.names = F, col.names = F)
  write.table(pvalsRange, paste(dataType,"_",mouseLine,"pvalsRange.txt",sep = ""), sep = " ",row.names = F, col.names = F)
  
}
pvals
# stopImplicitCluster()
pvalsRange
