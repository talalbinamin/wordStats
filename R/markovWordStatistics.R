#' Markov word statistics
#'
#' Calculate the word statistics (mean, variance, word count) using the Markov model.
#' @param seq The input DNA sequence
#' @param word_length The length of the words
#' @keywords word statistics, Bernoulli, motif finding
#' @return  returns the following in the form of a list \item{var}{the variance of the words} 
#' \item{mu}{the mean of the words}\item{N_w}{the word count}
#' @export calculateOverlapIndicator
#' @examples seq <- readDNAStringSet(filePath,format="fasta")
#'word_length <- 2
#'t <- markovWordStatistics(seq,word_length)
#'

markovWordStatistics <- function(seq,word_length,markovOrder) {
  
  markovModel <- markovModel(seq,markovOrder)
  noOfNtuples = 4^word_length
  myNames<-mkAllStrings(c("A","C","G","T"),word_length)
  onf<-oligonucleotideFrequency(seq,word_length+1,simplify.as="collapse")
  onfM1=rep(0,4^word_length)
  names(onfM1)<-myNames
  freqWord <-  oligonucleotideFrequency(seq, word_length)
  # find the mean and variance values of words for markov model
  mu = freqWord
  var = freqWord
  
  for(i in 1:noOfNtuples)
  {
   
    mu[i] <-meanMarkov(names(onfM1[i]),rep(0.25,4),markovModel$stationary,
                               markovModel$transition) 
    
    var[i] <- varianceMarkov(names(onfM1[i]),mu[i],width(seq),markovModel$stationary,
                                     markovModel$transition,exact=TRUE)
  
  }
  
  markov_word_statistics = list(mu=mu,var=var,N_w=freqWord)
  
  return(markov_word_statistics)
    
  
}