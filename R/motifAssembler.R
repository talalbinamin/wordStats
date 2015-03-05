motifAssembler <-
function(kmerVector,startKmer,position=0,exactCore=FALSE,extend=2,kmerWeightVector=NULL)
{
	#Recursion stop
	if(abs(position)>extend)
	{
		return(NULL)
	}
	else
	{
		#length of kmer
		k<-nchar(startKmer)
		
		leftMotif=kmerVector[substr(kmerVector,2,k)==substr(startKmer,1,k-1)][1]
		rightMotif=kmerVector[substr(kmerVector,1,k-1)==substr(startKmer,2,k)][1]
		
		#Each kmer that occurs in the motif is removed from the list
		if(position==0)
		{
			kmerVector=kmerVector[-which(kmerVector==startKmer)]
		}
		
		if(is.na(leftMotif)||(exactCore==TRUE&&position<0))
		{
			leftMotif=NULL
		}
		else		#Each kmer that occurs in the motif is removed from the list
		{
			kmerVector=kmerVector[-which(kmerVector==leftMotif)]
		}
		
		if(is.na(rightMotif)||(exactCore==TRUE&&position>0))
		{
			rightMotif=NULL
		}
		else		#Each kmer that occurs in the motif is removed from the list
		{
			kmerVector=kmerVector[-which(kmerVector==rightMotif)]
		}
		
		#Das Motif wird gespeichert als (kmer,position) pro Zeile
		
		#Special case, position==0, startKmer has to be stored in the motif
		if(position==0)
		{
			motifExtTMP<-c(startKmer,leftMotif,rightMotif)
			motifPosTMP<-c(rep(position,1),rep(position+1,length(leftMotif)),rep(position-1,length(rightMotif)))
		}
		else
		{
			motifExtTMP<-c(leftMotif,rightMotif)
			motifPosTMP<-c(rep(position+1,length(leftMotif)),rep(position-1,length(rightMotif)))
		}
		my_motif<-rbind(motifExtTMP,motifPosTMP,kmerWeightVector[motifExtTMP])
		if(position>=0&&!is.null(leftMotif))
		{
			my_motif<-cbind(my_motif,motifAssembler(kmerVector,leftMotif,position=position+1,exactCore=exactCore,extend=extend,kmerWeightVector=kmerWeightVector))
		}
		
		if(position<=0&&!is.null(rightMotif))
		{
			my_motif<-cbind(my_motif,motifAssembler(kmerVector,rightMotif,position=position-1,exactCore=exactCore,extend=extend,kmerWeightVector=kmerWeightVector))
		}
		return(my_motif)
	}
}
