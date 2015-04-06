#' Bernoulli word statistics
#'
#' Calculate the word statistics (mean, variance, word count) using the Bernoulli model.
#' @param seq The input DNA sequence
#' @param wordLength The length of the words
#' @keywords word statistics, Bernoulli, motif finding
#' @return  returns the following in the form of a list \item{var}{the variance of the words} 
#' \item{mu}{the mean of the words}\item{N_w}{the word count}
#' @export calculateOverlapIndicator 
#' @examples seq <- readDNAStringSet(filePath,format="fasta")
#'wordLength <- 2
#'t <- bernoulliWordStatistics(seq,wordLength)
#'

bernoulliWordStatistics <- function(seq,wordLength) {
  
  # a faster version for calculating the BErnoullis background model
  # 11/2/2015
  # tba
  
  lengthSeq = width(seq)
  numberOfLetters = 4 # A,C,G,T
  noOfNtuples = numberOfLetters^wordLength
  
  # compute the word counts (N) 
  
  freqLetter <- oligonucleotideFrequency(seq, 1,simplify.as="collapsed")
  freqWord <- oligonucleotideFrequency(seq, wordLength)
  
  probLetter = freqLetter/sum(freqLetter)
  
  wordNames <- colnames(freqWord)   
  
  
  # calculate the overlap factor (epsilon)
  
  epsilon <- sapply(wordNames, FUN = function(x) {
    return(unlist(calculateOverlapIndicator(x,x)))
  })
  
  epsilon <- t(epsilon)  
  
  wordVarianceOverlap <- freqWord
  probWord <- freqWord
  
  letterNames <- colnames(probLetter)
  
  names(freqLetter)=c("A","C","G","T")
  
  for (kk in 1:noOfNtuples){
    
    word = wordNames[kk]
    wordSplit=strsplit(word,"")[[1]]
    probWord[kk] = prod(probLetter[wordSplit[1:wordLength]])     
    term1 = (lengthSeq - wordLength + 1) * probWord[kk]*(1-probWord[kk]) 
    
    term2 = 0
    
    for (dd in 1: (wordLength - 1))
      
    {
      
      noOfTerms = (lengthSeq - wordLength - (dd-1))
      letters = substr(word,wordLength-abs(dd)+1,wordLength)
      wordSplit=strsplit(word,"")[[1]]
      prodTerm = prod(probLetter[wordSplit[(wordLength-dd+1):wordLength]])      
      
      term2Cov <- noOfTerms*(epsilon[kk,wordLength-dd] * prodTerm - probWord[kk])
      
      term2 = term2 + term2Cov 
      
    }
    
    term2 = term2  * 2*probWord[kk]
    
    wordVarianceOverlap[kk] = term1 + term2
    
  }
  
  return(list(var = wordVarianceOverlap, mu=probWord,word_count = freqWord))
}