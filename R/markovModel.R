markovModel<-function(seq,mo)
{
  myNames<-mkAllStrings(c("A","C","G","T"),mo)
  transition<-matrix(rep(0,(4^mo)^2),ncol=4^mo,nrow=4^mo)
  colnames(transition)=myNames
  rownames(transition)=myNames
  onf<-oligonucleotideFrequency(seq,mo+1,simplify.as="collapse")
  onf_m1=rep(0,4^mo)
  names(onf_m1)<-myNames
  for(i in 1:length(onf))
  {
    suffix<-substr(names(onf[i]),1,mo)
    onf_m1[suffix]=onf_m1[suffix]+onf[i]
  }
  
  for(oligo in 1:length(onf))
  {
    i<-substr(names(onf[oligo]),1,mo)
    j<-substr(names(onf[oligo]),2,2+mo-1)
    transition[i,j]<-onf[oligo]/onf_m1[i]
  }
  
  stationary<-(transition%*%transition%*%transition%*%transition%*%transition%*%transition%*%transition%*%transition%*%transition%*%transition%*%transition%*%transition%*%transition%*%transition)[1,]
  model=list(stationary=stationary,transition=transition,onf_m1=onf_m1)
  return(model)
}
