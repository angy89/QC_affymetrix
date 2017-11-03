# Required libraries
library(affy)
library(affyPLM)
library(simpleaffy)

#setwd(fileFolder)
print("Reading cel files...")
celFiles = list.files(fileFolder) #List of the names of cel files to be preprocessed
celFiles = paste(fileFolder,celFiles,sep="")
nFiles = length(celFiles) #number of cel files to analyse

# The RLE, NUSE and RNA degradation indices will be evaluated on random partition of the pool of samples

if(nFiles > maxFiles){
  groupSize = trunc(nFiles/maxFiles) #how many groups of files will be randomly created
  ElemToSplit = trunc(length(celFiles)/groupSize) * groupSize
  RS = split(celFiles[1:ElemToSplit],sample(1:groupSize))
  
  if(groupSize>1){
    RS[[groupSize]] = c(RS[[groupSize]],celFiles[(ElemToSplit+1):length(celFiles)])
  }
  
  save(RS,file=paste(RDATAfileFolder,"randomSubsample.RData",sep=""))
}


if(isParallel){
  
  library(foreach)
  library(doMC)
  registerDoMC(nCores)
  
  cat("FitPLM\n")
  #length(RS)
  PSetList = foreach(i=1:length(RS))%dopar% {
    message(i)
    affyB = ReadAffy(filenames = RS[[i]])
    Pset <- fitPLM(affyB)
  }
  
  save(PSetList,file=paste(RDATAfileFolder,"affy_pset_random_lists.RData",sep=""))
  
  
  cat("RNA DEg\n")
  RNADeg = foreach(i=1:length(RS),.combine = "cbind")%dopar% {
    affyB = ReadAffy(filenames = RS[[i]])
    # Pset <- fitPLM(affyB)
    deg = AffyRNAdeg(affyB)
    png(paste(RNADeg_fig_folder,"RNADeg_batch",i,".png",sep=""))
    plotAffyRNAdeg(deg)
    dev.off()
    summaryAffyRNAdeg(deg)
  }
  
  cat("RLE stat\n")
  RLE_stat = foreach(i = 1:length(PSetList),.combine = "cbind")%dopar% {
    Pset = PSetList[[i]]
    png(paste(RLE_fig_folder,"RLE_batch",i,".png",sep=""))
    rleP = RLE(Pset)
    dev.off()
    RLE(Pset,type="stat")
  }
  
  cat("NUSE stat\n")
  NUSE_stat = foreach(i = 1:length(PSetList),.combine = "cbind")%dopar% {
    Pset = PSetList[[i]]
    png(paste(NUSE_fig_folder,"NUSE_batch",i,".png",sep=""))
    nuseP = NUSE(Pset)
    dev.off()
    NUSE(Pset,type="stat")
  }
  
}else{
  RLE_stat = c()
  RLE_values = c()
  NUSE_stat = c()
  NUSE_values = c()
  RNADeg = c()
  
  print("Starting the analyses...")
  if(groupSize>1){
    pb = txtProgressBar(min = 1, max = length(RS),style=3)
  }
  
  PSetList = list()
  affyBList = list()
  for(i in 1:length(RS)){
    cat("Reading data...\n")
    affyB = ReadAffy(filenames = RS[[i]])
    cat("AffyPLM...\n")
    
    Pset <- fitPLM(affyB)
    PSetList[[i]] = Pset
    affyBList[[i]] = affyB
    
    #   # Relative Log Expression (RLE) values. Specifically,
    #   #these RLE values are computed for each probeset by comparing the expression value
    #   #on each array against the median expression value for that probeset across all arrays.
    cat("RLE...\n")
    
    png(paste(RLE_fig_folder,"RLE_batch",i,".png",sep=""))
    rleP = RLE(Pset)
    dev.off()
    RLE_stat = cbind(RLE_stat,RLE(Pset,type="stat"))
    # RLE_values = cbind(RLE_values,RLE(Pset,type="value"))
    
    #   #Normalized Unscaled Standard Errors (NUSE) can also be used for assessing quality. In
    #   #this case, the standard error estimates obtained for each gene on each array from fitPLM
    #   #are taken and standardized across arrays so that the median standard error for that
    #   #genes is 1 across all arrays
    cat("NUSE...\n")
    
    png(paste(NUSE_fig_folder,"NUSE_batch",i,".png",sep=""))
    nuseP = NUSE(Pset)
    dev.off()
    NUSE_stat = cbind(NUSE_stat,NUSE(Pset,type="stat"))
    #NUSE_values = cbind(NUSE_values,NUSE(Pset,type="value"))
    
    cat("DEG...\n")
    
    #measure the difference in the signal at the 5' and 3' of a gene
    deg = AffyRNAdeg(affyB)
    degS = summaryAffyRNAdeg(deg)
    png(paste(RNADeg_fig_folder,"RNADeg_batch",i,".png",sep=""))
    plotAffyRNAdeg(deg)
    dev.off()
    
    RNADeg = cbind(RNADeg,degS)
    if(groupSize>1){
      setTxtProgressBar(pb,i)
    }
  }
  
  if(groupSize>1){
    close(pb)
  }
  save(PSetList,affyBList,file=paste(RDATAfileFolder,"affy_pset_random_lists.RData",sep=""))
  
}


save(RNADeg,RLE_stat,NUSE_stat,file=paste(RDATAfileFolder,"outliers_stats.RData",sep=""))

png(paste(RNADeg_fig_folder,"RNADeg_slope.png",sep=""))
boxplot(RNADeg["slope",])
dev.off()

png(paste(RLE_fig_folder,"RLE_median.png",sep=""))
boxplot(RLE_stat["median",])
dev.off()

png(paste(NUSE_fig_folder,"NUSE_median.png",sep=""))
boxplot(NUSE_stat["median",])
dev.off()

print("Computing outliers...")

RLEout <- names(boxplot.stats(RLE_stat["median",])$out)  # outlier values.
NUSEout <- names(boxplot.stats(NUSE_stat["median",])$out)  # outlier values.
deg_bstat = boxplot(RNADeg["slope",])$stats[4] #value of the 4% quantile
DEGout = colnames(RNADeg)[which(RNADeg["slope",]>deg_bstat)]
outs = union(DEGout,union(RLEout,NUSEout)) #names of all the outliers samples

#M is a matrix of samples that contains a 1 or a 0 if the sample is considered outliers for one of the 6 measures
M = matrix(0,nrow=length(outs),ncol=3)
rownames(M) = outs
colnames(M) = c("RLE","NUSE","DEG")
M[outs %in% RLEout,"RLE"] = 1
M[outs %in% NUSEout,"NUSE"] = 1
M[outs %in% DEGout,"DEG"] = 1

M = cbind(M,rowSums(M))
colnames(M)[4] = "SUM"

outs_for_at_least_one_method = outs
outs = rownames(M)[M[,4]>1]
save(M,outs,outs_for_at_least_one_method,file=paste(RDATAfileFolder,"outliers.RData",sep=""))

print("End...")


