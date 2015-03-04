#' Calculate overlap between two words
#'
#' This function calculates the overlap between two words
#' @param 
#' @keywords 
#' @export
#' @examples
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
