% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/markovWordStatistics.R
\name{markovWordStatistics}
\alias{markovWordStatistics}
\title{Markov word statistics}
\usage{
markovWordStatistics(seq, word_length, markov_order)
}
\arguments{
\item{seq}{The input DNA sequence}

\item{word_length}{The length of the words}
}
\value{
returns the following in the form of a list \item{var}{the variance of the words}
\item{mu}{the mean of the words}\item{N_w}{the word count}
}
\description{
Calculate the word statistics (mean, variance, word count) using the Markov model.
}
\examples{
seq <- readDNAStringSet(filePath,format="fasta")
word_length <- 2
t <- markovWordStatistics(seq,word_length)
}
\keyword{Bernoulli,}
\keyword{finding}
\keyword{motif}
\keyword{statistics,}
\keyword{word}

