#these functions calculate all kinds of different statistics for detecting differential genes

calmoststatistics<-function(sample,n1,n2){
  result<-array(0,dim=c(dim(sample)[1],n2))
  temp<-generateorder(n2)
  a<-temp[,1]; b<-temp[,2];
  n<-n1+n2; samplex<-sample[,1:n1]; sampley<-sample[,(n1+1):n]
  medx<-apply(samplex,1,median)
  medy<-apply(sampley,1,median)
  diffxmedx<-abs(samplex-medx); diffymedy<-abs(sampley-medy)
  stdforeachrow<-apply(cbind(diffxmedx,diffymedy),1,median)+0.01
  for(ktest in 1:n2){
     if (ktest>1){meanoutlier<-apply(sampley[,1:ktest],1,mean)}
	else {meanoutlier<-sampley[,1]}
     #result[,ktest]<-((meanoutlier-medx)*ktest-a[ktest])/stdforeachrow #whether -a?
     result[,ktest]<-((meanoutlier-medx)*ktest/stdforeachrow/1.4826-a[ktest])/b[ktest]
  }
  mystatistics<-apply(result,1,max)
}

callssvstatistics<-function(sample,n1,n2)
  {
  result<-array(0,dim=c(dim(sample)[1],n2))
  varall<-array(0,dim=c(dim(sample)[1],n2))
  mystatistics<-array(0,dim=c(dim(sample)[1]))
  n<-n1+n2; samplex<-sample[,1:n1]; sampley1<-sample[,(n1+1):n]
  sampley<-t(apply(sampley1,1,sort,decreasing=T))
  meanx<-apply(samplex,1,mean)
  varx<-apply((samplex-meanx)^2,1,sum)
  for(ktest in 1:n2){
     if (ktest>1 && ktest<(n2-1))
     {
     meanoutlier<-apply(sampley[,1:ktest],1,mean)
     varoutlier<-apply((sampley[,1:ktest]-meanoutlier)^2,1,sum)
     meanother<-apply(sampley[,(ktest+1):n2],1,mean)
     varother<-apply((sampley[,(ktest+1):n2]-meanother)^2,1,sum)
     varall[,ktest]<-varoutlier+varother
     }
     if(ktest==1) 
     {
     meanoutlier<-sampley[,1]
     meanother<-apply(sampley[,(ktest+1):n2],1,mean)
     varother<-apply((sampley[,(ktest+1):n2]-meanother)^2,1,sum)
     varall[,ktest]<-varother
     }
     if(ktest==(n2-1))
     {
     meanoutlier<-apply(sampley[,1:ktest],1,mean)
     varoutlier<-apply((sampley[,1:ktest]-meanoutlier)^2,1,sum)
     varall[,ktest]<-varoutlier
     }
     if(ktest==n2)
     {
     meanoutlier<-apply(sampley[,1:ktest],1,mean)
     varoutlier<-apply((sampley[,1:ktest]-meanoutlier)^2,1,sum)
     varall[,ktest]<-varoutlier
     }
     result[,ktest]<-((ktest)*(meanoutlier-meanx)/(sqrt((varx+varall[,ktest]))+0.01))
     #result[,ktest]<-((meanoutlier-meanx)/(sqrt((varx+varall[,ktest]))+0.01))
  }
  for(i in 1:dim(sample)[1])
  {
  temp<-varall[i,]
  tag<-which.min(temp)
  mystatistics[i]<-result[i,tag]
  } 
  mystatistics
}

calortstatistics<-function(sample,n1,n2){
  n<-n1+n2; samplex<-sample[,1:n1]; sampley<-sample[,(n1+1):n]
  threshold<-apply(samplex,MARGIN=1,FUN=quantile,probs=0.75)+apply(samplex,1,IQR)
  outlierflag<-(sampley>threshold)
  medx<-apply(samplex,1,median)
  medy<-apply(sampley,1,median)
  diffxmedx<-abs(samplex-medx); diffymedy<-abs(sampley-medy)
  stdforeachrow<-apply(cbind(diffxmedx,diffymedy),1,median)+0.01
  ortstatistics<-(apply(sampley*outlierflag,1,sum)-apply(outlierflag,1,sum)*medx)/stdforeachrow
}

calosstatistics<-function(sample,n1,n2){
  n<-n1+n2; samplex<-sample[,1:n1]; sampley<-sample[,(n1+1):n]
  threshold<-apply(sample,MARGIN=1,FUN=quantile,probs=0.75)+apply(sample,1,IQR)
  outlierflag<-(sampley>threshold)
  medx<-apply(samplex,1,median)
  medy<-apply(sampley,1,median)
  med<-apply(sample,1,median)
  osstdforeachrow<-apply(abs(sample-apply(sample,1,median)),1,median)+0.01
  osstatistics<-(apply(sampley*outlierflag,1,sum)-apply(outlierflag,1,sum)*med)/osstdforeachrow
}

calcopastatistics<-function(sample,n1,n2){
  n<-n1+n2; samplex<-sample[,1:n1]; sampley<-sample[,(n1+1):n]
  med<-apply(sample,1,median)
  q90<-apply(sampley,1,quantile,probs=0.9)
  copastdforeachrow<-apply(abs(sample-med),1,median)+0.01
  capastatistics<-(q90-med)/copastdforeachrow
}

caltstatistics<-function(sample,n1,n2){
  n<-n1+n2; samplex<-sample[,1:n1]; sampley<-sample[,(n1+1):n]
  meanx<-apply(samplex,1,mean); meany<-apply(sampley,1,mean);
  tvar<-apply((samplex-meanx)^2,1,sum)+apply((sampley-meany)^2,1,sum)
  tstatistics<-(meany-meanx)/(sqrt(tvar)+0.01)
}

