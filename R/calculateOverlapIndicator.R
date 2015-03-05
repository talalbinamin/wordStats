#' Calculate word overlap
#'
#' Given two words, calculate overlap between them by shifting the second word over the first one, one letter at a time over 
#' @param word1 The first word
#' @param word2 The second word
#' @keywords word, overlap
#' @return  \item {Outputs the overlap between the two words. The size of the overlap is the same as the input word with the longest length}
#' @examples word1 = 'ACGA'
#' word2 = 'TCGC'
#' overlap <-calculateOverlapIndicator(word1,word2)
#'


calculateOverlapIndicator<-function(word1,word2)
{
  k=nchar(word1)
  eps=rep(0,k)
  for(u in 1:k)
  {
    if(substr(word1,k-u+1,k)==substr(word2,1,u))
    {
      eps[u]=1
    }
  }
  return(eps)
}
