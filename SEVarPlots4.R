# Current code, uncomment for parameters
# !/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 6) {
  stop("Needing 6 arguments")
}

dataSource = myargs[1]
mouseLine = myargs[2]
distr = myargs[3]
param1=as.numeric(myargs[4])
param2=as.numeric(myargs[5])
plotRows=as.numeric(myargs[6])
# 
# 
# 
# install.heteroplasmy = T
# 
# # install kimura package if required
# if(install.heteroplasmy) {
#   install.packages("devtools")
#   library("devtools")
#   devtools::install_github("kostasgian21/heteroplasmy")
# }



library("devtools")
library(kimura)
library(heteroplasmy)
library(ggplot2)
# library(ggpubr)



set.seed(74)


# #initial parameters
# dataSource="real" # real synthetic
# mouseLine="LE" # HB LE Freyer Broz path
# distr="rnorm"#rbeta rlnorm kimura rnorm
# param1=0.5
# param2=0.1
# plotRows=1
seqStart=6 # min sample size for synthetic data
seqEnd=50 # max sample size for synthetic data
seqStep=2 # increase step for  sample size for synthetic data
nrep=10000 #number of bootstrap trials

if (dataSource=="real") {
  numSamples=1
}else{
  numSamples=20 # number of repetitions, used when synthetic data are tested
}




if (dataSource=="real") {
  if (mouseLine=="Freyer") {
    mouseData=mousedataFreyer
  }else if(mouseLine=="Broz"){
    mouseData=mousedataBroz
  }else if(mouseLine=="LE"){
    mouseData=mousedataLE
    # mouseData=mouseData[,sample(1:length(mouseData), 12, replace=F)]
  }else if(mouseLine=="HB"){
    mouseData=mousedataHB
  }else{
    mouseData <- read.table(mouseLine, sep = "\t",
                            header = T)
  }
  ns=seq(1, length(mouseData), by=1)
  distr="NA"
  x_axis_label="Sample instance"
}else{
  mouseLine=distr
  x_axis_label="Sample size"
}


#loop throughout different values of n
if (dataSource=="real") {
  ns=seq(1, ncol(mouseData), by=1)
  lengt=length(mouseData)
}else{
  ns=seq(seqStart, seqEnd, by=seqStep)
  lengt=length(ns)
}

resultMeans= data.frame(samplesVarNorm=numeric(),
                        samplesVar=numeric(),
                        samplesVar2=numeric(),
                        bootVar=numeric(),
                        mlVar=numeric(),
                        anSEVNormal=numeric(),
                        anSEVMom=numeric(),
                        anSEVKimura=numeric(),
                        bootSEV=numeric(),
                        MLSEV=numeric())


for (z in 1:length(ns)) {
  
  
  result= data.frame(samplesVarNorm=numeric(),
                     samplesVar=numeric(),
                     samplesVar2=numeric(),
                     bootVar=numeric(),
                     mlVar=numeric(),
                     anSEVNormal=numeric(),
                     anSEVMom=numeric(),
                     anSEVKimura=numeric(),
                     bootSEV=numeric(),
                     MLSEV=numeric())
  

  #numSamples repetitions for each value of n to obtain average's
  j=1
  while (j<=numSamples) {
    
    #sample data/heteroplasmies from a known distribution
    
    if (dataSource=="real") {
      X.1=mouseData[which(!is.na(mouseData[,z])),z]
    }else{
      n=ns[z]
    }
    
    
    if (distr=="rnorm" && dataSource!="real") {
      X.1 = rnorm(n,mean = param1,sd = param2) 
    }
    # 
    
    if (distr=="rbeta"&& dataSource!="real") {
      X.1 = rbeta(n,shape1 = param1,shape2 = param2)
    }
    
    if (distr=="rlnorm"&& dataSource!="real") {
      X.1=rlnorm(n, meanlog = param1, sdlog = param2)
    }
    
    # pars=test_kimura(X.1)
    if (distr=="kimura"&& dataSource!="real") {
      X.1=rkimura(n,p = param1,b = param1)
    }
    
    X.1[which(X.1>1)]=1
    X.1[which(X.1<0)]=0
    
    # MoM estimation of p and b
    pars=estimate_parameters(X.1)
    p=pars[[1]]
    b=pars[[2]]
    
  
    n=length(X.1)
    
    #mean and variance of the sample heteroplasmy
    sampleVar=var(X.1)
    
    normAnalyt=analyticVar(X.1,normal = T)
    hstatAnalyt=analyticVar(X.1,method = "hstatistic")
    kimuraAnalyt=analyticVarKimura(X.1,p,b)
    SEB= bootstrapVar(X.1)[1]
    
    # # bootstrapp mean var point estimator
    bootVarMean=bootstrapVar(X.1)[2]
    
    # calculated, but not shown
    jackVariance=jackVar(X.1)[1]
    jackMean=jackVar(X.1)[2]
    
    # ML estimation of p and b
    fit = estimate_parameters_ml(X.1)
    h0=fit[1]
    b2=fit[2]

    mlKimura=analyticVarKimura(X.1,h0,b2)
    
    # maxlik var point estimator
    mlKimuraMean=sum((X.1-h0)^2)/(length(X.1)-1)
    
    # result is the dataframe that holds the result for the i-th repetition (same n)
    result[nrow(result) + 1,] = c(sampleVar,sampleVar,sampleVar,bootVarMean,
                                  mlKimuraMean,normAnalyt,hstatAnalyt,kimuraAnalyt,SEB,mlKimura)

    j=j+1
    
  }
  
  
  # take the means for the numSamples repetitions for the same n and report
  # them in the resultMeans dataframe that hold the overall output of the comparisons
  resultMeans[nrow(resultMeans) + 1,] = c(mean(result[,1]),mean(result[,2]),
                                          mean(result[,3]),mean(result[,4]),
                                          mean(result[,5]),mean(result[,6]),
                                          mean(result[,7]),mean(result[,8]),
                                          mean(result[,9]),mean(result[,10]))
}

# ###
# extra=2
# ns=seq(seqStart, seqEnd+seqStep*extra, by=seqStep)
# for (i in 1:extra) {
#   resultMeans[nrow(resultMeans) + 1,] = c(0,0,0,0,0,0,0,0,0,0)
# }
# ###


# Plotting and styling
row.names(resultMeans) <- ns
cl  <- c("#D55E00", "#56B4E9", "#009E73", "#CC79A7","#0072B2")
names(resultMeans)[1] <- "Normal appr"
names(resultMeans)[2] <- "Analytic"
names(resultMeans)[3] <- "FittedKimura"
names(resultMeans)[4]<-"Boots"
names(resultMeans)[5]<-"Maxlik"

resultMeansB=resultMeans
nsB=ns
prints=list()

resultMeans=resultMeansB[order(resultMeansB$`Normal appr`),]
ns=nsB[1:nrow(resultMeans)]
resultMeans=resultMeans[1:nrow(resultMeans),]


cl=c("Normal appr",
     "Analytic (H-statistic)",
     "FittedKimura (MoM params)",
     "Bootstrapping",
     "FittedKimura (maxlik params)")

rownames(resultMeans)=ns

df=data.frame(method=character(),sample=character(),mean=numeric(),sde=numeric())
# xs=c(4,4.1,8,8.1,11,11.1,14,14.1,29,30,44)
xs=ns
cols=vector()


for (i in 1:nrow(resultMeans)) {
  df[nrow(df)+1,]=c("1Norm. approx.",ns[i],resultMeans[i,1],2*resultMeans[i,6])
  # cols=c(cols,'#999999')
  cols=c(cols,cl[1])
}
for (i in 1:nrow(resultMeans)) {
  df[nrow(df)+1,]=c("2Analytic (MoM)",ns[i],resultMeans[i,2],2*resultMeans[i,7])
  # cols=c(cols,'#E69F00')
  cols=c(cols,cl[2])
}
for (i in 1:nrow(resultMeans)) {
  df[nrow(df)+1,]=c("3Fitted Kimura",ns[i],resultMeans[i,3],2*resultMeans[i,8])
  # cols=c(cols,'#E62E00')
  cols=c(cols,cl[3])
}
for (i in 1:nrow(resultMeans)) {
  df[nrow(df)+1,]=c("4Bootstrapping",ns[i],resultMeans[i,4],2*resultMeans[i,9])
  # cols=c(cols,'#E62E00')
  cols=c(cols,cl[4])
}
for (i in 1:nrow(resultMeans)) {
  df[nrow(df)+1,]=c("5Maxlik",ns[i],resultMeans[i,5],2*resultMeans[i,10])
  # cols=c(cols,'#E62E00')
  cols=c(cols,cl[5])
}
df$sample <- as.character(df$sample)
df$sample <- factor(df$sample, levels=unique(df$sample))

# df$dose=as.factor(df$dose)
df$mean=as.numeric(df$mean)
df$sde=as.numeric(df$sde)
df2=df


cl2  <- c("#D55E00", "#56B4E9", "#009E73", "#CC79A7","#0072B2")

p4<- ggplot(df2, aes(x=sample, y=mean, group=method)) + 
  geom_errorbar(aes(ymin=mean-sde, ymax=mean+sde,color=cols),
                size=1.0, width=0.9,
                position=position_dodge(0.8)) +
  geom_point(position = position_dodge(0.8),color="darkgray") +
  # expand_limits(y=0) + 
  xlab("Sample instance") + 
  ylab("Variance")+
  # labs(title=paste(dataSource, "heteroplasmy data from",mouseLine, "mice"), x=x_axis_label, y = "Variance")+
  # labs(title=paste("heteroplasmy data from","Freyer et al."), x=x_axis_label, y = "Variance")+
  # labs(title=paste(dataSource, "heteroplasmy data from a","normal", "distribution"), x=x_axis_label, y = "Variance")+
  labs(title=paste("stdError(Var(h)) of", "heteroplasmy data from",mouseLine, "mice"), x=x_axis_label, y = "Variance")+
  scale_colour_discrete("Methods to calculate stdError(Var(h))",breaks=cl)+
  # theme(text = element_text(size = 20))+
  # theme(axis.text.x = element_text(size = 5))+
  theme_classic(base_size = 15)+
  theme(axis.text.x = element_text(size = 10))
  # theme(axis.text = element_text(size = 15))




png(paste(dataSource,distr,"_",mouseLine,".png",sep = ""), width=1300, height=450)
p4 + facet_wrap(~ sample, nrow =plotRows,scales = "free_x",strip.position="bottom")+
  theme(strip.text.x = element_blank())
dev.off()


