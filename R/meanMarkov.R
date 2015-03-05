meanMarkov<-function(word,freq,stationary,transition)
{
  if(is.null(stationary))
  {
    #return(0.25^nchar(word))
    return(prod(freq[strsplit(word,"")[[1]]]))
    
    #return(prod(freq^oligonucleotideFrequency(DNAString(word),1,as.prob=F)))
  }
  else
  {
    w1<-substr(word,1,1)
    p=stationary[w1]
    for(i in 1:(nchar(word)-1))
    {
      w1<-substr(word,i,i)
      w2<-substr(word,i+1,i+1)
      p=p*transition[w1,w2]
    }
    return(p)
  }
}