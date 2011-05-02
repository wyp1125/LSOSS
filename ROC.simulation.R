source("morefuns.R")
source("calstatisticsnew.R")
#change the following line to simulate under different parameters
n1<-20;n2<-20;mu<-2;k<-10;p<-1000
#generate data from null distribution, 0 means 0 lines is from alternative, i.e. all null
sample0<-generateexpression(n1,n2,p,mu,k,0)
#generate data from alternative distribution
sample1<-generateexpression(n1,n2,p,mu,k,p)

#most statistics
mystat0<-calmoststatistics(sample0,n1,n2)
mystat1<-calmoststatistics(sample1,n1,n2)

#by varying the thresholds, we can generate the ROC curve
threshold<-quantile(c(mystat0,mystat1),probs=c(seq(0.01,0.99,by=0.05)))

fpr<-c(); tpr<-c()
for (t in threshold){
  fpr<-c(fpr,sum(mystat0>=t)/p)
  tpr<-c(tpr,sum(mystat1>=t)/p)
}
par(cex=1.4)
plot(fpr,tpr,type="p",pch=6,cex=1,main="n=m=20, u=2, k=10",xlab="False positive rate",ylab="True positive rate")

#ortstatistics
ortstat0<-calortstatistics(sample0,n1,n2)
ortstat1<-calortstatistics(sample1,n1,n2)

threshold<-quantile(c(ortstat0,ortstat1),probs=c(seq(0.01,0.99,by=0.05)))

fpr<-c(); tpr<-c()
for (t in threshold){
  fpr<-c(fpr,sum(ortstat0>=t)/p)
  tpr<-c(tpr,sum(ortstat1>=t)/p)
}
lines(fpr,tpr,type="p",lty=3,cex=1,pch=2)

#os statistics
osstat0<-calosstatistics(sample0,n1,n2)
osstat1<-calosstatistics(sample1,n1,n2)

threshold<-quantile(c(osstat0,osstat1),probs=c(seq(0.01,0.99,by=0.05)))

fpr<-c(); tpr<-c()
for (t in threshold){
  fpr<-c(fpr,sum(osstat0>=t)/p)
  tpr<-c(tpr,sum(osstat1>=t)/p)
}
lines(fpr,tpr,type="p",lty=3,cex=1,pch=3)

#copa statistics
copastat0<-calcopastatistics(sample0,n1,n2)
copastat1<-calcopastatistics(sample1,n1,n2)

threshold<-quantile(c(copastat0,copastat1),probs=c(seq(0.01,0.99,by=0.05)))

fpr<-c(); tpr<-c()
for (t in threshold){
  fpr<-c(fpr,sum(copastat0>=t)/p)
  tpr<-c(tpr,sum(copastat1>=t)/p)
}
lines(fpr,tpr,type="p",lty=3,cex=1,pch=4)



#t statistics
tstat0<-caltstatistics(sample0,n1,n2)
tstat1<-caltstatistics(sample1,n1,n2)

threshold<-quantile(c(tstat0,tstat1),probs=seq(0.01,0.99,by=0.05))
fpr<-c(); tpr<-c()
for (t in threshold){
  fpr<-c(fpr,sum(tstat0>=t)/p)
  tpr<-c(tpr,sum(tstat1>=t)/p)
}
lines(fpr,tpr,type="p",lty=3,cex=1,pch=5)

#lssv statistics
lssvstat0<-callssvstatistics(sample0,n1,n2)
lssvstat1<-callssvstatistics(sample1,n1,n2)

threshold<-quantile(c(lssvstat0,lssvstat1),probs=seq(0.01,0.99,by=0.05))
fpr<-c(); tpr<-c()
for (t in threshold){
  fpr<-c(fpr,sum(lssvstat0>=t)/p)
  tpr<-c(tpr,sum(lssvstat1>=t)/p)
}
lines(fpr,tpr,type="p",lty=3,cex=1,pch=1)


legend(0.75,0.57,legend=c("lsoss","most","ort","os","copa","t"),pch=c(1,6,2,3,4,5),cex=0.9)




