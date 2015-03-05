#' Bernoulli word statistics
#'
#' Calculate the word statistics (mean, variance, word count) using the Bernoulli model.
#' @param seq The input DNA sequence
#' @param word_length The length of the words
#' @keywords word statistics, Bernoulli, motif finding
#' @return  returns the following in the form of a list \item{var}{the variance of the words} 
#' \item{mu}{the mean of the words}\item{N_w}{the word count}
#' @export calculateOverlapIndicator
#' @examples library(Biostrings)
#'seq <- readDNAStringSet(filePath,format="fasta")
#'word_length <- 2
#'t <- bernoulliWordStatistics(seq,word_length)
#'

bernoulliWordStatistics <- function(seq,word_length) {
  
  # a faster version for calculating the BErnoullis background model
  # 11/2/2015
  # tba
  
  #   # Read the DNA sequence 
  #   seq=readDNAStringSet(file=paste("C:/users/amint/Desktop/Seq Data/",
  #                                   "Visel_mmForebrain_p300_top2000_",
  #                                   "maxMask30_500bp_concatenated.fasta", sep=""),format="fasta")
  
  #   seq=readDNAStringSet(file=paste("/mnt/projects/amint/synbio/",
  #                                   "Visel_mmForebrain_p300_top2000_",
  #                                   "maxMask30_500bp_concatenated.fasta", sep=""),format="fasta")
  
  length_seq = width(seq)
  #   word_length <- 3
  no_letters = 4 # A,C,G,T
  no_ntuples = no_letters^word_length
  
  #  indexing individual elements in a DNA string   
  #  seq[[1]][1]
  
  # compute the word counts (N) 
  
  freq_letter <- oligonucleotideFrequency(seq, 1,simplify.as="collapsed")
  freq_word <- oligonucleotideFrequency(seq, word_length)
  
  #   # compute the letter and word probability (these are the ground truth)
    prob_letter = freq_letter/sum(freq_letter)
  #   prob_word = freq_word/sum(freq_word)
  
  word_names <- colnames(freq_word)
  
  
  #   # compute the word probability (this is an estimate from the background model)
  #   
  #   # for Bernoulli model
  #   
  #   # find Y(w)
  #   
  #   
  #   matched_list <- lapply(word_names, FUN = function(x) {
  #     return(vmatchPattern(x, seq))
  #   })
  #   
  #   # returns the start indices of the word in the sequence
  #   word_start <- lapply(matched_list, FUN = function(x) {
  #     return(unlist(start(x)))
  #   })
  #   # assign 1's to locations (indices) where the word exists
  #   Y <- matrix(0, no_ntuples,length_seq)
  #   
  #   for (ii in 1:no_ntuples ){    
  #     
  #     Y[ii, unlist(word_start[[ii]])] <- 1 
  #     
  #   }
  #   
  #   # mean and variance for Bernoulli model (eq 3.2 and 3.3)
  #   word_mean = rowMeans(Y, na.rm = FALSE, dims = 1) 
  #   word_variance = word_mean*(1-word_mean)
  #   
  #   # word counts from Bernoulli model (Same as freq_word)
  #   N_w = rowSums(Y, na.rm = FALSE, dims = 1) 
  
  # calculate the overlap factor (epsilon)
  
  epsilon <- sapply(word_names, FUN = function(x) {
    return(unlist(calculateOverlapIndicator(x,x)))
  })
  
  epsilon <- t(epsilon)
  
#   # expected counts from Bernoulli model for overlap dependance
#   
#   w <- (length_seq - word_length + 1) * prob_word
#   
  # variance of expected counts from Bernoulli model for overlap dependance
  
  word_variance_overlap <- freq_word
  prob_word <- freq_word
  #   word_variance_overlap2  <- rep(NA, no_ntuples)
  # term 1 (where i=j) appears l-k+1 times (diagonals of the matrix of i=1,2,..l-k+1 and j-1,2,...,l-k+1)
  
  
  #   # the number of times term 2 occurs (to statisfy abs(i-j) <= k-1)
  #   no_term2 = 2*((length_seq-word_length):(length_seq-word_length-(word_length-1)))
  
  letter_names <- colnames(prob_letter)
  
  names(freq_letter)=c("A","C","G","T")
  
  for (kk in 1:no_ntuples){
    
    word = word_names[kk]
    wordSplit=strsplit(word,"")[[1]]
    prob_word[kk] = prod(prob_letter[wordSplit[1:word_length]])     
    term_1 = (length_seq - word_length + 1) * prob_word[kk]*(1-prob_word[kk])
    #         term2 = 0
    #         
    #         for (dd in 1: (word_length - 1))
    #           
    #         {
    #           
    #           no_terms = 2*(length_seq - word_length - (dd-1))
    #           letters = substr(word,word_length-abs(dd)+1,word_length)
    #           wordSplit=strsplit(word,"")[[1]]
    #           prod_term = prod(freq_letter[wordSplit[(word_length-dd+1):word_length]])      
    #           
    #           term2_cov <- prob_word[kk]*epsilon[kk,word_length-dd] * prod_term - prob_word[kk]^2
    #           
    #           term2 = term2 + (term2_cov * no_terms)
    #           
    #         }
    #         
    #         word_variance_overlap[kk] = term_1[kk] + term2
    #         
    #       }
    #       
    
    term2 = 0
    
    for (dd in 1: (word_length - 1))
      
    {
      
      no_terms = (length_seq - word_length - (dd-1))
      letters = substr(word,word_length-abs(dd)+1,word_length)
      wordSplit=strsplit(word,"")[[1]]
      prod_term = prod(prob_letter[wordSplit[(word_length-dd+1):word_length]])      
      
      term2_cov <- no_terms*(epsilon[kk,word_length-dd] * prod_term - prob_word[kk])
      
      term2 = term2 + term2_cov 
      
    }
    
    term2 = term2  * 2*prob_word[kk]
    
    word_variance_overlap[kk] = term_1 + term2
    
    #     
    #     word_variance_overlap2[kk] <-wordVarIID(word,c("A","C","G","T"),length_seq,c(prob_letter[1],prob_letter[2],
    #                                                                                  prob_letter[3],prob_letter[4]),
    #                                             prob_word[kk])
    #     
  }
  
  return(list(var = word_variance_overlap, mu=prob_word,word_count = freq_word))
}