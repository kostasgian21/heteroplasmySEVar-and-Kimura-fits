# Current code, uncomment for parameters
# !/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 3) {
  stop("Needing 3 arguments")
}


numberOfMCMCsamples=as.numeric(myargs[1])
nSamplesNoisedData=as.numeric(myargs[2])
reps=as.numeric(myargs[3])

library(kimura)
library(heteroplasmy)
library(foreach)

# setwd("C:/Users/kgi062/Desktop/mydata/code/heteroplasmy")


# # comment bellow for not-parallel code (and also line 151). Chabge also %dopar% to %do% in line 124
# library(doParallel)
# library(parallel)
# cores=2
# cl=makeCluster(cores)
# 
# registerDoParallel(cl)
# 
# clusterEvalQ(cl, expr={library(kimura)
#   library(heteroplasmy)})





genBinnedData = function() {

  
  wonList=list()
  
  md2=vector()#fig 1
  md2=c(rep(0,8),rep(5,35),rep(15,20),rep(25,10),rep(35,5),rep(45,4))
  md2=md2/100
  
  md2=c(rep(0,8),runif(35,min=1,max=10),runif(20,min=10,max=20),runif(10,min=20,max=30),runif(5,min=30,max=40),runif(4,min=40,max=50))/100
  
  wonList[[length(wonList)+1]]=md2
  
  md2=vector()#fig 2A
  md2=c(rep(0,9),rep(0.005,3),rep(0.015,4),rep(0.025,3),rep(0.035,1),rep(0.045,3),rep(0.055,2),rep(0.065,2),
        rep(0.075,3),rep(0.085,2),rep(0.095,1),rep(0.105,4),rep(0.125,1),rep(0.145,1),rep(0.155,1),rep(0.185,1))
  
  md2=c(rep(0,9),runif(3,min=0,max=1),runif(4,min=1,max=2),runif(3,min=2,max=3),runif(1,min=3,max=4),runif(3,min=4,max=5),runif(2,min=5,max=6),runif(2,min=6,max=7),runif(3,min=7,max=8),runif(2,min=8,max=9),runif(1,min=9,max=10),runif(4,min=10,max=11),runif(1,min=12,max=13),runif(1,min=14,max=15),runif(1,min=15,max=16),runif(1,min=18,max=19))/100
  
  wonList[[length(wonList)+1]]=md2
  
  md2=vector()#fig 2C
  md2=c(rep(0,10),rep(0.005,6),rep(0.015,1),rep(0.025,3),rep(0.035,1),rep(0.045,3),rep(0.055,1),rep(0.065,1),
        rep(0.145,1),rep(0.155,1),rep(0.195,1),rep(0.245,1),rep(0.195,1),rep(0.325,1))
  
  md2=c(rep(0,10),runif(6,min=0,max=1),runif(1,min=1,max=2),runif(3,min=2,max=3),runif(1,min=3,max=4),runif(3,min=4,max=5),runif(1,min=5,max=6),runif(1,min=6,max=7),runif(1,min=14,max=15),runif(1,min=15,max=16),runif(1,min=19,max=20),runif(1,min=32,max=33))/100
  
  wonList[[length(wonList)+1]]=md2
  
  md2=vector()#fig 3A
  md2=c(rep(0,6),rep(1,7),rep(3,3),rep(5,6),rep(7,3),rep(9,3),
        rep(11,4),rep(13,3),rep(17,4),rep(19,2),rep(21,2))
  md2=md2/100
  
  md2=c(rep(0,6),runif(7,min=0,max=2),runif(3,min=2,max=4),runif(6,min=4,max=6),runif(3,min=6,max=8),runif(3,min=8,max=10),runif(4,min=10,max=12),runif(3,min=12,max=14),runif(4,min=16,max=18),runif(2,min=18,max=20),runif(2,min=20,max=22))/100
  
  wonList[[length(wonList)+1]]=md2
  
  md2=vector()#fig 3C
  md2=c(rep(0,7),rep(1,1),rep(3,1),rep(5,3),rep(7,3),
        rep(11,3),rep(13,2),rep(15,1),rep(17,1),rep(19,1),rep(21,1),rep(23,1),rep(51,1))
  md2=md2/100
  
  md2=c(rep(0,7),runif(1,min=0,max=2),runif(1,min=2,max=4),runif(3,min=4,max=6),runif(3,min=6,max=8),runif(3,min=10,max=12),runif(2,min=12,max=14),runif(1,min=14,max=16),runif(1,min=16,max=18),runif(1,min=18,max=20),runif(1,min=20,max=22),runif(1,min=22,max=24),runif(1,min=50,max=52))/100
  
  wonList[[length(wonList)+1]]=md2
  
  md2=vector()#fig 4A
  md2=c(rep(0,24),rep(0.25,5),rep(0.75,5),rep(1.25,2),rep(1.75,6),rep(2.25,1),rep(2.75,3),
        rep(4.25,1),rep(5.25,2))
  md2=md2/100
  
  md2=c(rep(0,24),runif(5,min=0,max=0.5),runif(5,min=0.5,max=1),runif(2,min=1,max=1.5),runif(6,min=1.5,max=2),runif(1,min=2,max=2.5),runif(3,min=2.5,max=3),runif(1,min=4,max=4.5),runif(2,min=5,max=5.5))/100
  
  wonList[[length(wonList)+1]]=md2
  
  
  md2=vector()#fig 4B
  md2=c(rep(0,16),rep(0.75,16),rep(1.25,4),rep(1.75,5),rep(2.25,2),rep(2.75,2),
        rep(3.75,1),rep(4.75,1),rep(5.25,1),rep(7.75,1))
  md2=md2/100
  
  md2=c(rep(0,16),runif(16,min=0.5,max=1),runif(4,min=1,max=1.5),runif(5,min=1.5,max=2),runif(2,min=2,max=2.5),runif(2,min=2.5,max=3),runif(1,min=3.5,max=4),runif(1,min=4.5,max=5),runif(2,min=5,max=5.5),runif(1,min=7.5,max=8))/100
  
  wonList[[length(wonList)+1]]=md2
  
  
  md2=vector()#fig 5A
  md2=c(rep(0,4),rep(0.25,3),rep(1.25,1),rep(1.75,6),rep(2.25,3),rep(2.75,2),rep(3.25,3),
        rep(4.25,2),rep(4.75,1),rep(5.75,2),rep(6.25,2),rep(11.75,1),rep(12.25,1))
  md2=md2/100
  
  md2=c(rep(0,4),runif(3,min=0,max=0.5),runif(1,min=1,max=1.5),runif(6,min=1.5,max=2),runif(3,min=2,max=2.5),runif(2,min=2.5,max=3),runif(3,min=3,max=3.5),runif(2,min=4,max=4.5),runif(1,min=4.5,max=5),runif(2,min=5.5,max=6),runif(2,min=6,max=6.5),runif(1,min=11.5,max=12),runif(1,min=12.5,max=13))/100
  
  wonList[[length(wonList)+1]]=md2
  
  
  md2=vector()#fig 5C
  md2=c(rep(0,13),rep(0.75,4),rep(1.25,3),rep(1.75,1),rep(2.25,5),rep(2.75,5),rep(3.25,5),
        rep(3.75,1),rep(4.25,1),rep(4.75,1),rep(5.25,1),rep(5.75,2),rep(7.25,1),rep(8.75,1),rep(9.25,2))
  md2=md2/100
  
  md2=c(rep(0,13),runif(4,min=0.5,max=1),runif(3,min=1,max=1.5),runif(1,min=1.5,max=2),runif(5,min=2,max=2.5),runif(5,min=2.5,max=3),runif(5,min=3,max=3.5),runif(1,min=3.5,max=4),runif(1,min=4,max=4.5),runif(1,min=4.5,max=5),runif(1,min=5,max=5.5),runif(2,min=5.5,max=6),runif(1,min=7,max=7.5),runif(1,min=8.5,max=9),runif(2,min=9,max=9.5))/100
  
  wonList[[length(wonList)+1]]=md2
  
  
  md2=vector()#fig 6A*
  md2=c(rep(1,14),rep(6,9),rep(15,8),rep(25,4),rep(35,1),rep(45,3),
        rep(75,1),rep(85,1),rep(94,3),rep(99,16))
  md2=md2/100
  
  md2=c(runif(14,min=0,max=2),runif(9,min=2,max=10),runif(8,min=10,max=20),runif(4,min=20,max=30),runif(1,min=30,max=40),runif(3,min=40,max=50),runif(1,min=70,max=80),runif(1,min=80,max=90),runif(3,min=90,max=98),runif(16,min=98,max=100))/100
  
  wonList[[length(wonList)+1]]=md2
  
  md2=vector()#fig 6C
  md2=c(rep(6,1),rep(15,6),rep(25,14),rep(35,21),rep(45,8),
        rep(55,4))
  md2=md2/100
  
  md2=c(runif(1,min=2,max=10),runif(6,min=10,max=20),runif(14,min=20,max=30),runif(21,min=30,max=40),runif(8,min=40,max=50),runif(4,min=50,max=60))/100
  
  wonList[[length(wonList)+1]]=md2
  
  
  md2=vector()#fig 6E
  md2=c(rep(35,1),rep(45,11),rep(55,12),rep(65,17),rep(75,8),
        rep(85,1))
  md2=md2/100
  
  md2=c(runif(1,min=30,max=40),runif(11,min=40,max=50),runif(12,min=50,max=60),runif(17,min=60,max=70),runif(8,min=70,max=80),runif(1,min=80,max=90))/100
  
  wonList[[length(wonList)+1]]=md2
  
  md2=vector()#fig 7A
  md2=c(rep(1,3),rep(6,7),rep(15,27),rep(25,16),rep(35,4),
        rep(45,2))
  md2=md2/100
  
  md2=c(runif(3,min=0,max=2),runif(7,min=2,max=10),runif(27,min=10,max=20),runif(16,min=20,max=30),runif(4,min=30,max=40),runif(2,min=40,max=50))/100
  
  wonList[[length(wonList)+1]]=md2
  
  
  md2=vector()#fig 7C* 
  md2=c(rep(15,1),rep(25,1),rep(35,6),rep(45,10),rep(55,7),rep(65,5),rep(75,1))
  md2=md2/100
  
  md2=c(runif(1,min=10,max=20),runif(1,min=20,max=30),runif(6,min=30,max=40),runif(10,min=40,max=50),runif(7,min=50,max=60),runif(5,min=60,max=70),runif(1,min=70,max=80))/100
  
  wonList[[length(wonList)+1]]=md2
  
  md2=vector()#fig 7E* 
  md2=c(rep(45,3),rep(55,4),rep(65,5),rep(75,11),rep(85,9),rep(94,8),rep(99,12))
  md2=md2/100
  
  md2=c(runif(3,min=40,max=50),runif(4,min=50,max=60),runif(5,min=60,max=70),runif(11,min=70,max=80),runif(9,min=80,max=90),runif(8,min=90,max=98),runif(12,min=98,max=100))/100
  
  wonList[[length(wonList)+1]]=md2
  
  md2=vector()#fig 8A
  md2=c(rep(1,5),rep(4.5,7),rep(7.5,10),rep(12.5,7),rep(17.5,3),rep(22.5,5),
        rep(27.5,4),rep(32.5,2),rep(37.5,1))
  md2=md2/100
  
  md2=c(runif(5,min=0,max=2),runif(7,min=2,max=5),runif(10,min=5,max=10),runif(7,min=10,max=15),runif(3,min=15,max=20),runif(5,min=20,max=25),runif(4,min=25,max=30),runif(2,min=30,max=35),runif(1,min=35,max=40))/100
  
  wonList[[length(wonList)+1]]=md2

  return(wonList)
}






# #initial parameters
# numberOfMCMCsamples=50
# nSamplesNoisedData=2
reps=3


wonList=genBinnedData()

pvals= data.frame(Moments=numeric(),
                  maxlik=numeric())
pvalsRange=data.frame(Min=numeric(),
                      Max=numeric())

pvalsList=list()
pvalsRangeList=list()

for (k in 1:reps) {
  wonList=genBinnedData()
for (z in 1:16) {
  X.1=wonList[[z]]
  # X.1=abs(heteroplasmyShift(X.1,0.75))
  X.1[which(X.1>1)]=1
  X.1[which(X.1<0)]=0
  

  n=length(X.1)
  pvalstmp= data.frame(Moments=numeric(),
                       maxlik=numeric())
  mm <- foreach (i=1:nSamplesNoisedData, .combine=rbind) %do% {
    aa=kimura::test_kimura(X.1,num_MC = numberOfMCMCsamples)
    # bb=test_kimura_par(X.1,p=estimate_parameters(X.1)[1],b=estimate_parameters(X.1)[2],num_MC = numberOfMCMCsamples)
    bb=test_kimura(X.1,num_MC = numberOfMCMCsamples)
    c(aa$p.value,bb$p.value)
  }
  pvalstmp[1:nSamplesNoisedData,]=mm[1:nSamplesNoisedData,]
  pvals[z,]=c(mean(pvalstmp[,1]),mean(pvalstmp[,2]))
  
  # range of MoM pvals
  pvalsRange[z,]=c(min(pvalstmp[,1]),max(pvalstmp[,1]))

}
  pvalsList[[k]]=pvals
  # stopImplicitCluster()
  pvalsRangeList[[k]]=pvalsRange
}


pvals= data.frame(Moments=numeric(),
                  maxlik=numeric())
pvalsRange=data.frame(Min=numeric(),
                      Max=numeric())

for (i in 1:16) {
  s1=0
  s2=0
  minP=1.1
  maxP=-1
  for (j in 1:reps) {
    s1=s1+pvalsRangeList[[j]][[1]][i]
    s2=s2+pvalsRangeList[[j]][[2]][i]
    if (pvalsRangeList[[j]][[1]][i]<minP) {
      minP=pvalsRangeList[[j]][[1]][i]
    }
    if (pvalsRangeList[[j]][[1]][i]>maxP) {
      maxP=pvalsRangeList[[j]][[2]][i]
    }
  }
  pvals[i,]=c(s1/reps,s2/reps)
  pvalsRange[i,]=c(minP,maxP)
}
pvals
pvalsRange

# stopImplicitCluster()

write.table(pvals, "wonnapinijBinnedDataPvals.txt", sep = " ",row.names = F, col.names = F)
write.table(pvalsRange, "wonnapinijBinnedDataPRange.txt", sep = " ",row.names = F, col.names = F)


