## Individual

##########################
# Inidividual measurements
# Somewhere here I would have to filter dataNHim to retain only the ten with longest blades

for (col in c( "Year", "plotNo")) 
  dataNHim[,col] <- factor(dataNHim[,col])

  tail(dataNHim)
  dim(dataNHim)
  tail(dataNHpi)

for (col in c("bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter"))  
  dataNHim[,col] <- as.numeric(dataNHim[,col])

#### dataNHim order of plotNo needs to be the same as that in dataNHpi, They were sorted alphabetically before saving to the rdata object
################## MH add this
####  str(dataNHim)  

#### dataNHim2<-dataNHim[order(dataNHim$Order_In_Plot),]
################## the order of plot in msXim/dataNHim matter $$$$

# Complicated by the fact that there are no data for some im plots
emZpl <- model.matrix( ~ -1 + plotNo, data=dataNHim)
  dim(emZpl)
  #rownames(dataNHim)<-paste0("im_",1:nrow(dataNHim))
  #rownames(emZpl)<-rownames(dataNHim)
  #nlevels(dataNHim$plotNo)
  #nlevels(dataNHpi$plotNo)
  #dataNHpi$plotNo
  #dataNHim$plotNo
levels(dataNHim$plotNo) <- levels(dataNHpi$plotNo)
##OR USE THIS:  dataNHim$plotNo <- factor(dataNHim$plotNo, levels=levels(dataNHpi$plotNo))

msXim <- model.matrix( ~ -1 + plotNo, data=dataNHim)  # nrow=nrow(dataNHim), ncol=number of plotNo=nlevels(plotNo)
#rownames(msXim)<-rownames(dataNHim)
  dim(msXim)
  dim(dataNHim)
# row(msXim) is the individual measurements, ordered alphabeticaly in dataNHim; col(msXim) is the levels of plotNo matched to dataNHpi$plotNo
# row(msZ) is the plotsNo, ordered alphabetically in dataNHpi; col(msZ) is the order of fndr,FG,MG,plots(SP progeny) in vector u, and in hMat
# row(msXdb) is the plotsNo in dataNHpi,col(msXdb) are the different factors
# rownames(dataNHpi)<-dataNHpi$plotNo  
# msXdb <- model.matrix( ~ densityBlades + line + block + popChk, data=dataNHpi) ### adding densityBlades, msXdb
# 

msZim <- emZsp <- msXim %*% msZ
msXim <- msXim %*% msXdb
  dim(msXdb)
  dim(msXim)
  dim(msZim)
  dim(dataNHim)

msOutBLh <- mixed.solve(y=log(dataNHim$bladeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBMWh <- mixed.solve(y=log(dataNHim$bladeMaxWidth+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBTh <- mixed.solve(y=log(dataNHim$bladeThickness+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSLh <- mixed.solve(y=log(dataNHim$stipeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSDh <- mixed.solve(y=log(dataNHim$stipeDiameter+1), Z=msZim, K=hMat, X=msXim, SE=T)
h2hMatim <- c(heritability(msOutBLh), heritability(msOutBMWh), heritability(msOutBTh), heritability(msOutSLh), heritability(msOutSDh)) #
names(h2hMatim) <- c("bladeLength", "bladeMaxWidth",  "bladeThickness", "stipeLength", "stipeDiameter") 

  round(h2hMatim,3)

#2019Yr  
#bladeLength  bladeMaxWidth bladeThickness    stipeLength  stipeDiameter 
#0.008          0.106          0.384          0.274          0.464 

#2020Yr
# bladeLength  bladeMaxWidth bladeThickness    stipeLength  stipeDiameter 
#  0.175          0.079          0.766          0.555          0.820 
 
#BothYr  
  #bladeLength  bladeMaxWidth bladeThickness    stipeLength  stipeDiameter 
  #0.059          0.079          0.530          0.420          0.656 
  
  
names(h2hMat) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")
names(h2hMatim) <- c("bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter")

allBLUPs <- cbind(msOutWWPh$u, msOutDWPMh$u, msOutPDWh$u, msOutDBh$u, msOutBLh$u, msOutBMWh$u,  msOutBTh$u, msOutSLh$u, msOutSDh$u)
colnames(allBLUPs) <- c("WWP", "DWpM", "PDW", "BDns", "BLen", "BMax",  "BThk", "SLen", "SDia")
#rownames(allBLUPs) is the same as that in u, and hMat

pdf(paste0("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/","CorrPlot_",yr,".pdf"))
corrplot::corrplot.mixed(cor(allBLUPs), diag="n", tl.cex=0.9, tl.col=1)
dev.off()



##### Obtaining the BLUPs of the GPs to estimat its BV  
sampledSP <- which(nchar(rownames(allBLUPs)) < 11)  # The fndr SP is max of 10 characters
releaseGP <- nchar(rownames(allBLUPs))   
releaseGP <- which(10 < releaseGP & releaseGP < 15)  #FG and MG from fndr SP is max of 14 characters
progenySP <- which(nchar(rownames(allBLUPs)) > 15)   #the Cross name is larger than 15 characters

fndNames <- rownames(allBLUPs)[sampledSP]
hasDesc <- sapply(fndNames, function(n) grep(n, rownames(allBLUPs)[progenySP]))
nDesc <- sapply(hasDesc, length)

# Different colors for fndr, GP, progeny SP crosses
colSPGP <- c(rep("black", length(sampledSP)), rep("dark green", length(releaseGP)), rep("dark red", length(progenySP)))
colSPGP[which(nDesc == 0)] <- "grey"

pdf(paste0("DWpMvsIndividual",yr,".pdf"))
plot(allBLUPs[,2], pch=16, col=colSPGP, xlab="Individual number", ylab="Dry weight per plot") # !!The second col is DWpM
dev.off()

rnPos <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]>0))] #!! The second col is DWpM
rnNeg <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]<0))] #!! The second col is DWpM
locPos <- table(substring(rnPos, 6, 7))
locNeg <- table(substring(rnNeg, 6, 7))

allBLUPsDF <- as.data.frame(allBLUPs)
allBLUPsDF <- cbind(pedigree=rownames(allBLUPs), allBLUPsDF)
#rownames(allBLUPsDF)[progenySP] <- paste("SL",yr,"-UCONN-S", dataNHpi$plotNo[1:length(progenySP)], sep="")  #### Add the plot No to it
rownames(allBLUPsDF)[progenySP]<-dataNHpi$crossID[1:length(progenySP)]

pdf(paste0("DWpMvsPDWsporophytes",yr,".pdf"))
plot(allBLUPsDF$PDW[progenySP], allBLUPsDF$DWpM[progenySP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Sporophytes", cex.lab=1.3)
dev.off()
pdf(paste0("DWpMvsPDWgametophytes",yr,".pdf"))
plot(allBLUPsDF$PDW[releaseGP], allBLUPsDF$DWpM[releaseGP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Gametophytes", cex.lab=1.3)
dev.off()

# Add gametophyte sex as a column to the dataframe
getGPsex <- function(gpName){
  sexNum <- strsplit(gpName, "-")[[1]][4]
  sex <- substring(sexNum, 1, 1)
  return(ifelse(sex %in% c("F", "M"), sex, NA))
}
gpSex <- sapply(rownames(allBLUPsDF), getGPsex)  ### Apply to the list
gpSex

# A selection index for sporophytes and gametophytes that weights DWpM twice as high as PDW
allBLUPsDF <- cbind(pedigree=allBLUPsDF$pedigree, gpSex=gpSex, index=2*allBLUPsDF$DWpM + allBLUPsDF$PDW, allBLUPsDF[,-1]) #Reorder cols
saveRDS(allBLUPsDF, file=paste0("allBLUPsDF_analyzeNH_",yr,".rds"))

# Best sporophytes
nTop<-20  ### Pick the top 20
bestSP_DWpM <- allBLUPsDF[progenySP,][order(allBLUPsDF$DWpM[progenySP], decreasing=T)[1:nTop],] # Order based on DWpM
bestSP_DWpM[,-(1:2)] <- bestSP_DWpM[,-(1:2)] %>% round(3) # Only keep the numeric part
write.table(bestSP_DWpM, paste0("BestSPbyDryWgtPerM_",yr,".txt"), quote=F, row.names=T, col.names=T)

bestSP_idx <- allBLUPsDF[progenySP,][order(allBLUPsDF$index[progenySP], decreasing=T)[1:nTop],]
bestSP_idx[,-(1:2)] <- bestSP_idx[,-(1:2)] %>% round(3)
write.table(bestSP_idx, paste0("BestSPbyDWpMandPDW_",yr,".txt"), quote=F, row.names=T, col.names=T)

# Best gametophytes
bestGP_DWpM <- allBLUPsDF[releaseGP,][order(allBLUPsDF$DWpM[releaseGP], decreasing=T)[1:nTop],]
bestGP_DWpM[,-(1:2)] <- bestGP_DWpM[,-(1:2)] %>% round(3)
write.table(bestGP_DWpM[,-1], paste0("BestGPbyDryWgtPerM_",yr,".txt"), quote=F, row.names=T, col.names=T)
bestGP_idx <- allBLUPsDF[releaseGP,][order(allBLUPsDF$index[releaseGP], decreasing=T)[1:nTop],]
bestGP_idx[,-(1:2)] <- bestGP_idx[,-(1:2)] %>% round(3)
write.table(bestGP_idx[,-1], paste0("BestGPbyDWpMandPDW_",yr,".txt"), quote=F, row.names=T, col.names=T)

# Ordered female and male gametophytes
temp <- allBLUPsDF[allBLUPsDF$gpSex%in%"F",-(1:2)]  ### The F and M is levels, so cannot use "=="
temp <- temp[order(temp$DWpM, decreasing=T),] %>% round(3)
write.table(temp, paste0("FemaleGP_OrderedByDryWgtPerM_",yr,".txt"), quote=F, row.names=T, col.names=T)

temp <- allBLUPsDF[allBLUPsDF$gpSex%in%"M",-(1:2)]
temp <- temp[order(temp$DWpM, decreasing=T),] %>% round(3)
write.table(temp, paste0("MaleGP_OrderedByDryWgtPerM_",yr,".txt"), quote=F, row.names=T, col.names=T)


