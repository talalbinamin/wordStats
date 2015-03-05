consensusMatrixFromAlignment <-
function(motifAlignment,as.pwm=FALSE,unknownWeight=0.025,kmerWeightVector=NULL)
{
l=nchar(motifAlignment[1])
pfm=matrix(0,ncol=l,nrow=4)
rownames(pfm)=c("A","C","G","T")
if(is.null(kmerWeightVector))
{
for(i in 1:l)
{
letterPosI=sapply(motifAlignment,substr,i,i)
pfm["A",i]=sum(letterPosI=="A")+unknownWeight*sum(letterPosI=="-")
pfm["C",i]=sum(letterPosI=="C")+unknownWeight*sum(letterPosI=="-")
pfm["G",i]=sum(letterPosI=="G")+unknownWeight*sum(letterPosI=="-")
pfm["T",i]=sum(letterPosI=="T")+unknownWeight*sum(letterPosI=="-")

}
}
else
{
#weightsTMP=kmerWeightVector[gsub("-","",motifAlignment)]
for(i in 1:l)
{
letterPosI=sapply(motifAlignment,substr,i,i)
pfm["A",i]=sum(as.integer(letterPosI=="A")*kmerWeightVector)+unknownWeight*sum(as.integer(letterPosI=="-")*kmerWeightVector)
pfm["C",i]=sum(as.integer(letterPosI=="C")*kmerWeightVector)+unknownWeight*sum(as.integer(letterPosI=="-")*kmerWeightVector)
pfm["G",i]=sum(as.integer(letterPosI=="G")*kmerWeightVector)+unknownWeight*sum(as.integer(letterPosI=="-")*kmerWeightVector)
pfm["T",i]=sum(as.integer(letterPosI=="T")*kmerWeightVector)+unknownWeight*sum(as.integer(letterPosI=="-")*kmerWeightVector)

}

pfm[pfm<=0]=10^(-10)
}
if(as.pwm==TRUE)
{
pfm=t(t(pfm)/colSums(pfm))
}
return(pfm)
}
