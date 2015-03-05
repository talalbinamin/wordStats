printMotif <-
function(motif,file="",format="alignment",unknownWeight=0.025)
{
	
	motifStart<-min(as.integer(motif[2,]))
	motifEnd<-max(as.integer(motif[2,]))
	outMotif<-NULL
    kmerWeightVector=NULL
    if(nrow(motif)==3)
    {
        kmerWeightVector=as.double(motif[3,])
	}
    j<-1
	for(i in 1:ncol(motif))
	{
		outMotif[j]=paste(paste(rep("-",motifEnd-as.integer(motif[2,i])),collapse="",sep=""),motif[1,i],paste(rep("-",as.integer(motif[2,i])-motifStart),collapse="",sep=""),sep="")
		j<-j+1
	}
	if(file!="")
	{
		write(outMotif,sep="",file=file)
	}
	
	if(format=="PFM")
	{
		outMotif<-consensusMatrixFromAlignment(outMotif,kmerWeightVector=kmerWeightVector,unknownWeight=unknownWeight)	
	}
	else if(format=="PWM")
	{
		outMotif<-consensusMatrixFromAlignment(outMotif,as.pwm=T,kmerWeightVector=kmerWeightVector,unknownWeight=unknownWeight)	
	}
	else if(format=="consensus")
	{
		outMotif<-consensusMatrixFromAlignment(outMotif,kmerWeightVector=kmerWeightVector,unknownWeight=unknownWeight)	
		outMotif<-paste((rownames(outMotif)[apply(outMotif,2,which.max)]),sep="",collapse="")
	}
	
	return(outMotif)
}
