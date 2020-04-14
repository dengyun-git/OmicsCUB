#### In this short sample code, we aim to explore relationship between two variables, Pi (theoretical distribution density) and Freq (observed distribution density)

#### One gene subsequence configuration corresponds to one 'Pi' value and one 'Freq' value

#### main body
library("ggplot2")
source("procDta.R")
source("reduceDta.R")

###(1) read data frame from file 'species1.csv', 'aa,L,Pi,GeneID  
dta <-read.csv('species1.csv',stringsAsFactors=F,header=T,sep=",")

###(2) Each aa of each length has a fitting, then each aa has a file to record fitting results
tmp3 <- read.csv(text="aa,L,Pi,Freq") ### empty data frame

###call function 'procDta': process raw dataset 'dta' to 'tmp2' for fitting, 'tmp2' contains 'one type aa, one type L, Pi,Freq, GeneID'
for (aa in c('Ile','Gly','Arg')){
  
  tmp2 <- procDta(dta,aa,5) 
  
  ###(3) simple linear regression
  ft<-lm(Pi,Freq,data=tmp2)
  hist(ft.fit$resid,main=paste("Histogram of Residuals of", aa))     ### check residul distribution
  slp<-(summary(ft)$coefficients)[2]                                 ### slope of fitted line
  err<- summary(ft)$adj.r.squared                                   ###adjusted R squred error
  
  ###(4) write the fitting result to file containing fittings for all the lengths 
  write(file=paste('/Users/yd46/',aa,"Slope.csv",sep=''),c(5,slp,err),append="TRUE",sep=",")
  
  tmp3<-rbind(tmp3,tmp2) ### tmp3 combine values for all aa
}  
###(5) visulise relationship between 'Pi' and 'Freq' for different aa  
glt <- ggplot(tmp3,aes(x=Pi,y=Freq,group="aa")) + geom_point(size=2,aes(colour="aa")) + 
  geom_smooth(method='lm',se=F) + scale_colour_discrete( ) 
glt + labs(title="fitting", x="Pi", y="Freq")

################# functions called by main body####################################
#Calculate empirical frequencies for each unique Pi and Retruns a normalised subset ready for fitting
procDta<-function(dta,aa,leng){
  tmp2<-dta[which(dta[["aa"]]==aa & dta[["L"]]==leng),] #### tmp2 contain one type aa and one type length
  
  l3 <- matrix(rep(0, each = nrow(tmp2)),nrow=nrow(tmp2),ncol=1)
  l3 <- reduceDta(tmp2)
  colnames(l3)<-c("Freq")
  
  tmp2<-cbind(tmp2,l3)
  tmpSum<- sum(tmp2[,"Freq"])   ### all counts for one length
  tmp2[,"Freq"]<- tmp2[,"Freq"]/tmpSum  ### normalised empirical frequency
  tmp2
}

reduceDta<-function(sset){
  un<-unique(sset[["Pi"]])
  for (i in un){
    l3[which(sset[["Pi"]]==i),1]<-length(sset[which(sset[["Pi"]]==i),]) ### find empirical count for each unique configuration (namely unique Pi)
  }
  l3 
}
