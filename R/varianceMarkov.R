varianceMarkov<-function(word,mu_word,l,stationary,transition,exact=TRUE)
{
  k=nchar(word)
  eps=calculateOverlapIndicator(word,word)
#   mu_word=prob(word,stationary=stationary,transition=transition)
  wordSplit=strsplit(word,"")[[1]]
  #first term for covariance calculation
  term_1=(l-k+1)*mu_word*(1-mu_word)
  term_2=0
  for(d in 1:(k-1))
  {
    prod_2=1
    for(j in (k-d+1):k)
    {
      prod_2=prod_2*transition[wordSplit[j-1],wordSplit[j]]
    }
    term_2=term_2+(l-k-d+1)*(eps[k-d]*prod_2-mu_word)
  }
  term_2=2*mu_word*term_2
  
  #third term
  term_3=0
  if(exact==TRUE)
  {
    matrix_pow_t=transition
    for(t in 1:(l-2*k+1))
    {
      term_3=term_3+(l-2*k-t+2)*(1/(stationary[wordSplit[1]])*matrix_pow_t[wordSplit[k],wordSplit[1]]-1)
      matrix_pow_t_tmp=matrix_pow_t%*%transition
      if(identical(matrix_pow_t,matrix_pow_t_tmp))
      {break}
      matrix_pow_t=matrix_pow_t_tmp
      
    }
    term_3=term_3*mu_word*mu_word*2
  }
  
  return(term_1+term_2+term_3)
}