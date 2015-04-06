varianceMarkov<-function(word,muWord,l,stationary,transition,exact=TRUE)
{
  k=nchar(word)
  eps=calculateOverlapIndicator(word,word)
#   muWord=prob(word,stationary=stationary,transition=transition)
  wordSplit=strsplit(word,"")[[1]]
  #first term for covariance calculation
  term1=(l-k+1)*muWord*(1-muWord)
  term2=0
  for(d in 1:(k-1))
  {
    prod2=1
    for(j in (k-d+1):k)
    {
      prod2=prod2*transition[wordSplit[j-1],wordSplit[j]]
    }
    term2=term2+(l-k-d+1)*(eps[k-d]*prod2-muWord)
  }
  term2=2*muWord*term2
  
  #third term
  term3=0
  if(exact==TRUE)
  {
    matrixPowT=transition
    for(t in 1:(l-2*k+1))
    {
      term3=term3+(l-2*k-t+2)*(1/(stationary[wordSplit[1]])*matrixPowT[wordSplit[k],wordSplit[1]]-1)
      matrixPowTTmp=matrixPowT%*%transition
      if(identical(matrixPowT,matrixPowTTmp))
      {break}
      matrixPowT=matrixPowTTmp
      
    }
    term3=term3*muWord*muWord*2
  }
  
  return(term1+term2+term3)
}