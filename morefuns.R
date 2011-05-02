# This file contains two functions, generateorder use simulation to generate 
# the mean and standard deviation used in standardizing M_{ik} for differentk.
#The other function generateonesim is used to generate simulated gene 
#expression data, 

generateorder<-function(n){
  sample<-rnorm(1000*n)
  msample<-matrix(sample,nrow=1000)
  ordered<-apply(msample,MARGIN=1,FUN=sort,decreasing=T)
  cumsumordered<-apply(ordered,2,cumsum)
  a<-apply(cumsumordered,1,mean)
  cumsumordered<-cumsumordered-a
  b<-apply(cumsumordered,1,sd)
  cbind(a,b)
}

#simulate one set with p genes and n=n1+n2 arrays, mu is the magnitude of 
#difference for two samples, k is number of outliers, each row in the disease
#sample is finally sorted, h1 is number of rows that are DE
generateexpression<-function(n1,n2,p,mu,k,h1=1){
  n<-n1+n2 #number of arrays
  samplex<-matrix(rnorm(p*n1),nrow=p)
  sampley<-matrix(rnorm(p*n2),nrow=p)
  if (h1>=1) { sampley[1:h1,1:k]=sampley[1:h1,1:k]+mu }
  sampley<-t(apply(sampley,1,sort,decreasing=T))
  cbind(samplex,sampley)
}
