### Both Years
###################################
# Plot information analysis
###################################

#plot(dataNHpiBoth_C[dataNHpiBoth_C$Year==2019,]$densityBlades~dataNHpiBoth_C[dataNHpiBoth_C$Year==2019,]$dryWgtPerM)

library(rrBLUP)
rm(list=ls())
load("dataNHpi_withChk_3_sets_PhotoScore23.rdata")

load("dataNHim_withChk_3_sets_PhotoScore0123.rdata")

ls()
head(dataNHpiBoth_C)

#colMeans(dataNHpiBoth_C[dataNHpiBoth_C$PhotoScore>1,20:26],na.rm=TRUE) ### !!!! USE PHOTO SCORE 2 & 3 !

### !!!!
dataNHpi<-dataNHpiBoth_C  
dataNHim<-dataNHimboth_C
yr<-"Both"


load(paste0("hMat_PedNH_CCmat_fndrMrkData_",yr,"_PhotoScore23.rdata"))
dim(dataNHpi)
dim(hMat)
tail(dataNHpi)  
spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))    ### Progeny sps, col 2 and 3 are all values
nSp <- length(spRows)
nSp
dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "Z", substr(dataNHpi$plotNo, 1, 2), "ES")  # Checks VS ES

# Experimental design variables should be factors
for (col in c( "Year", "plotNo","GrowDays","femaPar", "femaParLoc", "malePar", "maleParLoc", "PlantingDens", "block", "line","popChk")) 
  dataNHpi[,col] <- factor(dataNHpi[,col])

nlevels(dataNHpi$plotNo)

#MHuang edit this:
# withinLoc: variable is 1 if the cross was between gametophytes sampled from
# the same location, zero otherwise
dataNHpi$withinLoc <- ifelse(as.vector(dataNHpi$femaParLoc) == as.vector(dataNHpi$maleParLoc), 1, 0) # WithinLoc is 1

dataNHpi$development[dataNHpi$development=="#N/A"] <-NA

# Make data cols numeric. Enforce data is numeric
for (col in c("wetWgtPlot", "lengthPlot", "wetWgtPerM","percDryWgt",  "dryWgtPerM","densityBlades")) 
  dataNHpi[,col] <- as.numeric(dataNHpi[,col])

# If percentDryWeigth is NA, then the plot WetWeightPerM and DryWeightperM should be set at NA, WetWeigthPerM may be just rope
dataNHpi[is.na(dataNHpi$percDryWgt), c("wetWgtPerM", "dryWgtPerM")] <- NA
keepRows <- !is.na(dataNHpi$percDryWgt)
keepRows

# isSelf: variable is 1 if the cross was between gametophytes from the same 
# sporophyte, zero otherwise
## MHuang edit
fndrF1<-strsplit(as.character(dataNHpi$femaPar), split="-", fixed=T)
fndrM1<-strsplit(as.character(dataNHpi$malePar),split="-",fixed=T)
fndrF<-sapply(fndrF1, function(vec) paste(vec[1:3], collapse="-"))
fndrM<-sapply(fndrM1, function(vec) paste(vec[1:3],collapse="-"))

isSelf<-ifelse(fndrF==fndrM,1,0)
dataNHpi$isSelf <- isSelf 
str(dataNHpi)

# Experimental design variables should be factors
for (col in c("Year","plotNo","GrowDays", "femaPar", "femaParLoc", "malePar", "maleParLoc", "line", "block", "date", "popChk", "withinLoc", "isSelf")) 
  dataNHpi[,col] <- as.factor(dataNHpi[,col])

### MHuang edit this
#dataNHpiES<-dataNHpi[dataNHpi$popChk=="ES",]
#dataNHpiChk<-dataNHpi[dataNHpi$popChk!="ES",]

#dataNHpi<-rbind(dataNHpiES,dataNHpiChk)  # 569 rows ES + 47 rows of Chk
dim(dataNHpi)
tail(dataNHpi)
### 1. Naive heritability using amatrix, Analyze GOM, BOTH years
msX <- model.matrix( ~ line%in%Year+block%in%Year+Year+popChk, data=dataNHpi)  ### !!!! nested effects???? DID NOT give the full Ranks of X
msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]

# Construction of Z assumes that all checks are at the bottom of the dataset
# Z has dimensions nrow(dataNHpi) x ncol(aMat)
# the first part all 0s, diag() identity matrix

# modify msZ to add the checks
# two same cross for a different plot would be having 

dim(msX)
colnames(msX)

### aMat has fndr+GP+ all the plot SPs in it.
### sort 
  dim(dataNHpi_RMchk)
  str(dataNHpi_RMchk)
  
dataNHpi_RMchk<-dataNHpi[!dataNHpi$crossID=="Check",] # Subset without chk levels
dataNHpi_RMchk$Crosses<-as.factor(as.character(dataNHpi_RMchk$Crosses)) # This RMed the check levels
  
nrowplot<-nrow(dataNHpi_RMchk)
nrowchk<-nrow(dataNHpi) - nrowplot

## 1. fndrs and GPs
fndr_GPMat<-matrix(0, nrowplot, nrow(aMat) - nSp)
  dim(fndr_GPMat)

## 2. plot level Sps
nSpMat<-model.matrix( ~ -1+Crosses,data=dataNHpi_RMchk)  # No intercept; sorted inside the 
  dim(nSpMat)
colnames(nSpMat)<-str_replace(colnames(nSpMat),"Crosses","")  

### the Crosses in nSpMat,make them the same order as those crosses in aMat
nSpMat<-nSpMat[,colnames(hMat[,(ncol(fndr_GPMat)+1):ncol(hMat)])]

## 3. chks rows  
chkSpMat<-matrix(0, nrowchk, nrow(aMat))
  dim(chkSpMat)  

## 4. make msZ
msZ <- rbind(cbind(fndr_GPMat,nSpMat),chkSpMat )
  dim(msZ)
  dim(aMat)
  
  which(msZ[19,]==1) #325
  which(msZ[129,]==1) #325
  hMat[325:326,325:326]

# Calculate heritability from the mixed.solve output
heritability <- function(msOut){
  return(msOut$Vu / (msOut$Vu + msOut$Ve))
}


########### 6. FINALE !! RUN Model WITHOUT blade density as covariate, but use hMat as the relationship matrx
write.csv(dataNHpi,paste0("dataNHpi_Last_Used_in_Model_",yr,".csv"))

msX <- model.matrix( ~ line%in%Year+block%in%Year+Year+popChk, data=dataNHpi)  ### !!!! nested effects???? DID NOT give the full Ranks of X
msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]  
msOutDBh <- mixed.solve(y=dataNHpi$densityBlades, Z=msZ, K=hMat, X=msX, SE=T)
msOutWWPh <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msX, SE=T)
msOutDWPMh <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msX, SE=T)
msOutPDWh <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msX, SE=T)

h2hMat <- c(heritability(msOutWWPh), heritability(msOutDWPMh), heritability(msOutPDWh),heritability(msOutDBh))
names(h2hMat) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt","densityBladesPerM")
  round(h2hMat,3) 

#    wetWgtPlot        dryWgtPerM        percDryWgt densityBladesPerM 
# 0.506             0.440             0.041             0.209  
 # wetWgtPlot        dryWgtPerM        percDryWgt densityBladesPerM 
  #0.502             0.437             0.044             0.205 

############# Individual
  for (col in c( "Year", "plotNo")) 
    dataNHim[,col] <- factor(dataNHim[,col])
  
  tail(dataNHim)
  dim(dataNHim)
  tail(dataNHpi)
  
  for (col in c("bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter"))  
    dataNHim[,col] <- as.numeric(dataNHim[,col])
  
  #### dataNHim of plotNo needs to be the same as that in dataNHpi, They were sorted alphabetically before saving to the rdata object
  ################## the order of plot in msXim/dataNHim matter $$$$
  
  # Complicated by the fact that there are no data for some im plots
  emZpl <- model.matrix( ~ -1 + plotNo, data=dataNHim)
  dim(emZpl)
  emZpl[1:5,1:6]
  #rownames(dataNHim)<-paste0("im_",1:nrow(dataNHim))
  #rownames(emZpl)<-rownames(dataNHim)
  #nlevels(dataNHim$plotNo)
  #nlevels(dataNHpi$plotNo)
  #dataNHpi$plotNo
  #dataNHim$plotNo
  #levels(dataNHim$plotNo) <- levels(dataNHpi$plotNo)
  ##OR USE THIS:  dataNHim$plotNo <- factor(dataNHim$plotNo, levels=levels(dataNHpi$plotNo))
  
  ########## Now the levels in dataNHim is more than those in dataNHpi, because the dataNHpi is only using plots with photo score 2,3. But dataNHim is having all plots that produced individual measurements
  
  dataNHim2<-subset(dataNHim,dataNHim$plotNo%in%dataNHpi$plotNo) #!!!!
  dim(dataNHim2)
  dim(dataNHim)
  dataNHim2$plotNo<-factor(dataNHim2$plotNo)
  nlevels(dataNHim2$plotNo) 
  nlevels(dataNHpi$plotNo)
  nlevels(dataNHim$plotNo)
  
  levels(dataNHim2$plotNo) <- levels(dataNHpi$plotNo)
  dataNHim<-dataNHim2
  
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
  #msX <- model.matrix( ~ line%in%Year+block%in%Year+Year+popChk, data=dataNHpi)
  msZim <- emZsp <- msXim %*% msZ
  #msXim <- msXim %*% msXdb
  msXim<-msXim%*%msX
  
  dim(msX)
  dim(msXim)
  dim(msZim)
  dim(dataNHim)
  msXim<-msXim[, apply(msXim, 2, function(v) !all(v == 0))]
  
  msOutBLh <- mixed.solve(y=log(dataNHim$bladeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
  msOutBMWh <- mixed.solve(y=log(dataNHim$bladeMaxWidth+1), Z=msZim, K=hMat, X=msXim, SE=T)
  msOutBTh <- mixed.solve(y=log(dataNHim$bladeThickness+1), Z=msZim, K=hMat, X=msXim, SE=T)
  msOutSLh <- mixed.solve(y=log(dataNHim$stipeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
  msOutSDh <- mixed.solve(y=log(dataNHim$stipeDiameter+1), Z=msZim, K=hMat, X=msXim, SE=T)
  h2hMatim <- c(heritability(msOutBLh), heritability(msOutBMWh), heritability(msOutBTh), heritability(msOutSLh), heritability(msOutSDh)) #
  names(h2hMatim) <- c("bladeLength", "bladeMaxWidth",  "bladeThickness", "stipeLength", "stipeDiameter") 
  
  round(h2hMatim,3)

  #BothYear_PhotoScore23
  #bladeLength  bladeMaxWidth bladeThickness    stipeLength  stipeDiameter 
  #0.060          0.115          0.480          0.219          0.564
  
  
  names(h2hMat) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")
  names(h2hMatim) <- c("bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter")
  
  allh2h<-c(h2hMat,h2hMatim)
  allh2h
  write.csv(allh2h,paste0("heritability_all_",yr,".csv"))
  
  allBLUPs <- cbind(msOutWWPh$u, msOutDWPMh$u, msOutPDWh$u, msOutDBh$u, msOutBLh$u, msOutBMWh$u,  msOutBTh$u, msOutSLh$u, msOutSDh$u)
  colnames(allBLUPs) <- c("WWP", "DWpM", "PDW", "BDns", "BLen", "BMax",  "BThk", "SLen", "SDia")
  #rownames(allBLUPs) is the same as that in u, and hMat
  
  pdf(paste0("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/","CorrPlot_",yr,"_PhotoScore23.pdf"))
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
  nDesc
  
  # Different colors for fndr, GP, progeny SP crosses
  colSPGP <- c(rep("black", length(sampledSP)), rep("dark green", length(releaseGP)), rep("dark red", length(progenySP)))
  colSPGP[which(nDesc == 0)] <- "grey"
  
  pdf(paste0("DWpMvsIndividual",yr,"_PhotoScore23.pdf"))
  plot(allBLUPs[,2], pch=16, col=colSPGP, xlab="Individual number", ylab="Dry weight per plot") # !!The second col is DWpM
  dev.off()
  
  rnPos <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]>0))] #!! The second col is DWpM
  rnNeg <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]<0))] #!! The second col is DWpM
  locPos <- table(substring(rnPos, 6, 7))
  locNeg <- table(substring(rnNeg, 6, 7))
  
  allBLUPsDF <- as.data.frame(allBLUPs)
  allBLUPsDF <- cbind(pedigree=rownames(allBLUPs), allBLUPsDF)
  #rownames(allBLUPsDF)[progenySP] <- paste("SL",yr,"-UCONN-S", dataNHpi$plotNo[1:length(progenySP)], sep="")  #### Add the plot No to it
  ### THIS DOES NOT WORK because it has "." rathe than "-", so later it is hard to strsplit()
  # rownames(allBLUPsDF)[progenySP]<-c(dataNHpi$crossID[1:length(progenySP)])
  head(allBLUPsDF)
  tail(allBLUPsDF)
  
  
  pdf(paste0("DWpMvsPDWsporophytes",yr,"_PhotoScore23.pdf"))
  plot(allBLUPsDF$PDW[progenySP], allBLUPsDF$DWpM[progenySP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Sporophytes", cex.lab=1.3)
  dev.off()
  pdf(paste0("DWpMvsPDWgametophytes",yr,"_PhotoScore23.pdf"))
  plot(allBLUPsDF$PDW[releaseGP], allBLUPsDF$DWpM[releaseGP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Gametophytes", cex.lab=1.3)
  dev.off()
  
  allBLUPsDF <- cbind(pedigree=allBLUPsDF$pedigree, index=2*allBLUPsDF$DWpM + allBLUPsDF$PDW, allBLUPsDF[,-1]) #Re
  cols
  saveRDS(allBLUPsDF, file=paste0("allBLUPsDF_analyzeNH_",yr,"_PhotoScore23.rds"))
  
  write.csv(allBLUPsDF,paste0("allBLUPs_DF_",yr,".csv"))
  
  
  # Best sporophytes
  nTop<-20  ### Pick the top 20
  bestSP_DWpM <- allBLUPsDF[progenySP,][order(allBLUPsDF$DWpM[progenySP], decreasing=T)[1:nTop],] # based on DWpM
  bestSP_DWpM[,-(1:2)] <- bestSP_DWpM[,-(1:2)] %>% round(3) # Only keep the numeric part
  write.table(bestSP_DWpM, paste0("BestSPbyDryWgtPerM_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)
  
  bestSP_idx <- allBLUPsDF[progenySP,][order(allBLUPsDF$index[progenySP], decreasing=T)[1:nTop],]
  bestSP_idx[,-(1:2)] <- bestSP_idx[,-(1:2)] %>% round(3)
  write.table(bestSP_idx, paste0("BestSPbyDWpMandPDW_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)
  
  # Best gametophytes
  bestGP_DWpM <- allBLUPsDF[releaseGP,][order(allBLUPsDF$DWpM[releaseGP], decreasing=T)[1:nTop],]
  bestGP_DWpM[,-(1:2)] <- bestGP_DWpM[,-(1:2)] %>% round(3)
  write.table(bestGP_DWpM[,-1], paste0("BestGPbyDryWgtPerM_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)
  bestGP_idx <- allBLUPsDF[releaseGP,][order(allBLUPsDF$index[releaseGP], decreasing=T)[1:nTop],]
  bestGP_idx[,-(1:2)] <- bestGP_idx[,-(1:2)] %>% round(3)
  write.table(bestGP_idx[,-1], paste0("BestGPbyDWpMandPDW_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)
  
  
