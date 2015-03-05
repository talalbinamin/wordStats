calculateMotifs <-
function(kmerWeightVector,searchSeq=NULL,max.test=5,motifNumber=3,extend=2,format="consensus",unknownWeight=0.01,plotCluster=FALSE,max.mismatch=1,...)
{
    kmerVector=names(sort(kmerWeightVector,d=T))
    motifSet<-list()
    for(i in 1:max.test)
    {
        startKmer<-kmerVector[i]
        motifSet[[i]]<-printMotif(motifAssembler(kmerVector,startKmer,position=0,exactCore=FALSE,extend=extend,kmerWeightVector=kmerWeightVector),file="",format="PWM",unknownWeight=unknownWeight)
        #startKmer=names(revComIndex[revComIndex[startKmer]])      
        #motifSet[i*2]<-printMotif(motifAssembler(kmerVector,startKmer,position=0,exactCore=FALSE),file="testRC.logo",format="consensus")
        names(motifSet)[i]=consensusString(motifSet[[i]],threshold=0.8,ambiguityMap="N")
        if(!is.null(searchSeq))
        {
            patternTMP=consensusString(motifSet[[i]],threshold=0.8,ambiguityMap="N")
            numNTMP=nchar(gsub("[^N]","",patternTMP,perl=T))
#tmp=matchPWM(pfmToPWM(motifSet[[i]]),searchSeq)
            tmp=matchPattern(patternTMP,searchSeq,max.mismatch=numNTMP+max.mismatch)
            if(length(tmp)>0)
            {
            pcm=consensusMatrix(tmp,baseOnly=T)[1:4,]
            motifSet[[i]]=t(t(pcm)/colSums(pcm))
            names(motifSet)[i]=consensusString(pcm,ambiguityMap="N")
            }
        }
    }

consensusMotifs=gsub("N","",names(motifSet),perl=T)
names(consensusMotifs)=names(motifSet)
    #Analysis of consensusmotifs:
    motifDist<-stringDist(consensusMotifs)

    motifClust<-hclust(motifDist)
    if(plotCluster==TRUE)
    {
        plot(motifClust,...)
    }
    motifClust$labels<-motifSet
    #plot(motifClust)
    return(motifSet[duplicated(cutree(motifClust,k=motifNumber))==FALSE])
}
