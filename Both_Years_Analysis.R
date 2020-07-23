### Both Years

###################################
# Plot information analysis
###################################
library(rrBLUP)

rm(list=ls())
load("dataNHpi_withChk_3_sets.rdata")
load("dataNHim_withChk_3_sets.rdata")


### !!!!
dataNHpi<-dataNHpiBoth_C  
dataNHim<-dataNHimboth_C
yr<-"Both"


load(paste0("hMat_PedNH_CCmat_fndrMrkData_",yr,".rdata"))
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
dataNHpiES<-dataNHpi[dataNHpi$popChk=="ES",]
dataNHpiChk<-dataNHpi[dataNHpi$popChk!="ES",]

dataNHpi<-rbind(dataNHpiES,dataNHpiChk)  # 569 rows ES + 47 rows of Chk
  dim(dataNHpi)
  tail(dataNHpi)
### 1. Naive heritability using amatrix, Analyze GOM, BOTH years
msX <- model.matrix( ~ line%in%Year+block%in%Year+Year+popChk, data=dataNHpi)  ### !!!! nested effects???? DID NOT give the full Ranks of X
msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]

# Construction of Z assumes that all checks are at the bottom of the dataset
# Z has dimensions nrow(dataNHpi) x ncol(aMat)
# the first part all 0s, diag() identity matrix

  dim(msX)
  colnames(msX)
msZ <- rbind(cbind(matrix(0, nSp, nrow(aMat) - nSp), diag(nSp)), matrix(0, nrow(dataNHpi) - nSp, nrow(aMat)))

# Calculate heritability from the mixed.solve output
heritability <- function(msOut){
  return(msOut$Vu / (msOut$Vu + msOut$Ve))
}

### Heritability for traits using pedigree-based relationship matrix
# h2naive means unconditional heritability

msOutWWP1 <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=aMat, X=msX)
msOutDB1 <- mixed.solve(y=dataNHpi$densityBlades, Z=msZ, K=aMat, X=msX)
msOutDWPM1 <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=aMat, X=msX)
msOutPDW1 <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=aMat, X=msX)
h2naive <- c(heritability(msOutDB1), heritability(msOutWWP1), heritability(msOutDWPM1), heritability(msOutPDW1))
names(h2naive) <- c("densityBlades", "wetWgtPlot", "dryWgtPerM", "percDryWgt")
  round(h2naive,3)
  
#densityBlades    wetWgtPlot    dryWgtPerM    percDryWgt 
 # 0.117         0.470         0.468         0.022  

  
  
# Heritability conditional on blade density using pedigree-based relationship
# Density Blades covariate
# The covariate relative to the focal trait could either remove sources of external random error or be a component that contributes genetic variation. The two are not mutually exclusive. To the extent that the former is prevalent the conditional heritability will be higher than the unconditional heritability. To the extent that the latter is prevalent, the reverse will hold. 

# 2. Add in the Density Blades covariate component, using A matrix:

msXdb <- model.matrix( ~ densityBlades + line%in%Year+block%in%Year+Year + popChk, data=dataNHpi) ### !!!! Change???
msXdb <- msXdb[, apply(msXdb, 2, function(v) !all(v == 0))]

msOutWWP <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=aMat, X=msXdb, SE=T)
msOutDWPM <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=aMat, X=msXdb, SE=T)
msOutPDW <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=aMat, X=msXdb, SE=T)

# extract the regression coefficient 
msOutDWPM$beta


h2covar <- c(heritability(msOutWWP), heritability(msOutDWPM), heritability(msOutPDW))
names(h2covar) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")

  round(h2covar,3)

  #wetWgtPlot dryWgtPerM percDryWgt 
  #0.368      0.376      0.000

#msXdb <- model.matrix( ~ densityBlades + line%in%Year+block%in%Year+Year + popChk, data=dataNHpi) ### !!!! Change???
### 3. Add in the Density Blades covariate component, BUT using H matrix:
msOutWWPh <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msXdb, SE=T)
msOutDWPMh <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msXdb, SE=T)
msOutPDWh <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msXdb, SE=T)

h2hMat <- c(heritability(msOutWWPh), heritability(msOutDWPMh), heritability(msOutPDWh))
names(h2hMat) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")
  round(h2hMat,3)
  
#wetWgtPlot dryWgtPerM percDryWgt 
#0.368      0.378      0.000 
  
### 4. Add in only the GrowDays covariate components, using H matrix:---- This is confounded with Year
# msX4 <- model.matrix( ~ GrowDays + line%in%Year+block%in%Year+Year + popChk, data=dataNHpi) ### !!!! Change???
# 
# msOutWWPh <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msX4, SE=T)
# msOutDWPMh <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msX4, SE=T)
# msOutPDWh <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msX4, SE=T)
# msOutDBh <- mixed.solve(y=dataNHpi$densityBlades, Z=msZ, K=hMat, X=msX, SE=T)
# h2hMat4 <- c(heritability(msOutWWPh), heritability(msOutDWPMh), heritability(msOutPDWh),heritability(msOutDBh))
# names(h2hMat4) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt","densityBlades")
#   round(h2hMat4,3)
# 

  ### 5. Using Development stage at painting as covariates, using H matrix
  
  msOutWWPh <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msX, SE=T)
  msOutDBh <- mixed.solve(y=dataNHpi$densityBlades, Z=msZ, K=hMat, X=msX, SE=T)
  msOutDWPMh <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msX, SE=T)
  msOutPDWh <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msX, SE=T)
  
  # The checks don't have this development measure, so assign at random and
  # do multiple times (nRepeats)
  
  nRepeats <- 2
  
  nIsNAdev <- sum(is.na(dataNHpi$development))
  fracDevelIs1 <- sum(dataNHpi$development == 1, na.rm=T) / sum(!is.na(dataNHpi$development))
  develBLUPs <- develVu <- develVe <- h2covarPSi <- NULL
  
  for (multImp in 1:nRepeats){
    dataNHpi$develImp <- dataNHpi$development
    dataNHpi$develImp[is.na(dataNHpi$development)] <- ifelse(runif(nIsNAdev) < fracDevelIs1, 1, 4)
    msOutDevelh <- mixed.solve(y=dataNHpi$develImp, Z=msZ, K=hMat, X=msX, SE=T)
    
    develBLUPs <- cbind(develBLUPs, msOutDevelh$u)
    develVu <- c(develVu, msOutDevelh$Vu)
    develVe <- c(develVe, msOutDevelh$Ve)
    
    # Conditional heritabilities: Paint SP devel
    msXps <- model.matrix( ~ develImp + line%in%Year + block%in%Year+Year + popChk, data=dataNHpi) ### !!!BothYears
    msXps<-msXps[, apply(msXps, 2, function(v) !all(v == 0))]
    msOutWWPps <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msXps, SE=T)
    
    msOutDWPMps <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msXps, SE=T)
    msOutPDWps <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msXps, SE=T)
    
    h2covarPS <- c(heritability(msOutWWPps), heritability(msOutDWPMps), heritability(msOutPDWps))
    names(h2covarPS) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")
    h2covarPSi <- cbind(h2covarPSi, h2covarPS)
    
  }
  
  msOutDevelh$u <- rowMeans(develBLUPs)
  msOutDevelh$Vu <- mean(develVu)
  msOutDevelh$Ve <- mean(develVe)
  h2covarPS <- rowMeans(h2covarPSi)
  
  h2hMat5 <- c(heritability(msOutWWPh), heritability(msOutDWPMh), heritability(msOutPDWh), heritability(msOutDBh), heritability(msOutDevelh))
  names(h2hMat5) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt", "densityBlades", "paintSPdevel")
  
  round(h2hMat5,3)
  
  #  wetWgtPlot    dryWgtPerM    percDryWgt densityBlades  paintSPdevel 
  #  0.369         0.357         0.014         0.219         0.512 
  
  # 2020 is not feasible, because the 2020 only has one level in development
  #
  
  # BothYr
  # wetWgtPlot    dryWgtPerM    percDryWgt densityBlades  paintSPdevel 
  # 0.470         0.468         0.008         0.134         0.257 
  
  
  
  
  
  
  
  
  
  
  
  
  
# 
# ##########################
# # Inidividual measurements
# # Somewhere here I would have to filter dataNHim to retain only the ten with longest blades
# 
# for (col in c( "Year", "plotNo")) 
#   dataNHim[,col] <- factor(dataNHim[,col])
# 
# tail(dataNHim)
# #CrossedSP<-read.csv("Crossed_WildSP_2019_2020.csv",sep=",",header=TRUE)
# dim(dataNHim)
# 
# for (col in c("bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter"))  
#   dataNHim[,col] <- as.numeric(dataNHim[,col])
# 
# #### !!!!! dataNHim order of plotNo needs to be the same as that in dataNHpi???
# ################## MH add this
# ####dataNHim$Order_In_Plot<-as.numeric(dataNHim$Order_In_Plot)
# ####  str(dataNHim)  
# 
# #### dataNHim2<-dataNHim[order(dataNHim$Order_In_Plot),]
# ################## Does the order of plot in msXim/dataNHim matter????????
# 
# # Complicated by the fact that there are no data for some im plots
# emZpl <- model.matrix( ~ -1 + plotNo, data=dataNHim)
# levels(dataNHim$plotNo) <- levels(dataNHpi$plotNo)
# dataNHim$plotNo <- factor(dataNHim$plotNo, levels=levels(dataNHpi$plotNo))
# 
# msXim <- model.matrix( ~ -1 + plotNo, data=dataNHim)  # nrow=nrow(dataNHim), ncol=number of plotNo=nlevels(plotNo)
# msZim <- emZsp <- msXim %*% msZ
# msXim <- msXim %*% msXdb
#   dim(msXdb)
#   dim(msXim)
#   dim(msZim)
#   dim(dataNHim)
# # bladeLength, bladeMaxWidth, bladeWid10cm, bladeThickness, stipeLength, stipeDiameter, stipeHollow, bullations, fouling
# 
# msOutBLh <- mixed.solve(y=log(dataNHim$bladeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
# msOutBMWh <- mixed.solve(y=log(dataNHim$bladeMaxWidth+1), Z=msZim, K=hMat, X=msXim, SE=T)
# msOutBTh <- mixed.solve(y=log(dataNHim$bladeThickness+1), Z=msZim, K=hMat, X=msXim, SE=T)
# msOutSLh <- mixed.solve(y=log(dataNHim$stipeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
# msOutSDh <- mixed.solve(y=log(dataNHim$stipeDiameter+1), Z=msZim, K=hMat, X=msXim, SE=T)
# h2hMatim <- c(heritability(msOutBLh), heritability(msOutBMWh), heritability(msOutBTh), heritability(msOutSLh), heritability(msOutSDh)) #
# names(h2hMatim) <- c("bladeLength", "bladeMaxWidth",  "bladeThickness", "stipeLength", "stipeDiameter") 
# 
#   round(h2hMatim,3)
# 
# 
# names(h2hMat4) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt", "densityBlades")
# names(h2hMatim) <- c("bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter")
# 
# allBLUPs <- cbind(msOutWWPh$u, msOutDWPMh$u, msOutPDWh$u, msOutDBh$u, msOutBLh$u, msOutBMWh$u,  msOutBTh$u, msOutSLh$u, msOutSDh$u)
# colnames(allBLUPs) <- c("WWP", "DWpM", "PDW", "BDns", "BLen", "BMax",  "BThk", "SLen", "SDia")
# 
# pdf(paste0("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/","CorrPlot_",yr,".pdf"))
# corrplot::corrplot.mixed(cor(allBLUPs), diag="n", tl.cex=0.9, tl.col=1)
# dev.off()
# 
# 
# 
# ##### Obtaining the BLUPs of the GPs to estimat its BV  
# sampledSP <- which(nchar(rownames(allBLUPs)) < 11)  #????
# releaseGP <- nchar(rownames(allBLUPs))   #????
# releaseGP <- which(10 < releaseGP & releaseGP < 15)  #????
# progenySP <- which(nchar(rownames(allBLUPs)) > 15)   #????
# 
# fndNames <- rownames(allBLUPs)[sampledSP]
# hasDesc <- sapply(fndNames, function(n) grep(n, rownames(allBLUPs)[progenySP]))
# nDesc <- sapply(hasDesc, length)
# colSPGP <- c(rep("black", length(sampledSP)), rep("dark green", length(releaseGP)), rep("dark red", length(progenySP)))
# colSPGP[which(nDesc == 0)] <- "grey"
# plot(allBLUPs[,2], pch=16, col=colSPGP, xlab="Individual number", ylab="Dry weight per plot")
# rnPos <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]>0))]
# rnNeg <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]<0))]
# locPos <- table(substring(rnPos, 6, 7))
# locNeg <- table(substring(rnNeg, 6, 7))
# 
# allBLUPsDF <- as.data.frame(allBLUPs)
# allBLUPsDF <- cbind(pedigree=rownames(allBLUPs), allBLUPsDF)
# rownames(allBLUPsDF)[progenySP] <- paste("SL18-UCONN-S", dataNHpi$plotNo[1:199], sep="") 
# 
# pdf("DWpMvsPDWsporophytes.pdf")
# plot(allBLUPsDF$PDW[progenySP], allBLUPsDF$DWpM[progenySP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Sporophytes", cex.lab=1.3)
# dev.off()
# pdf("DWpMvsPDWgametophytes.pdf")
# plot(allBLUPsDF$PDW[releaseGP], allBLUPsDF$DWpM[releaseGP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Gametophytes", cex.lab=1.3)
# dev.off()
# 
# # Add gametophyte sex as a column to the dataframe
# getGPsex <- function(gpName){
#   sexNum <- strsplit(gpName, "-")[[1]][4]
#   sex <- substring(sexNum, 1, 1)
#   return(ifelse(sex %in% c("F", "M"), sex, NA))
# }
# gpSex <- sapply(rownames(allBLUPsDF), getGPsex)
# gpSex
# # A selection index for sporophytes and gametophytes that weights DWpM twice as high as PDW
# allBLUPsDF <- cbind(pedigree=allBLUPsDF$pedigree, gpSex=gpSex, index=2*allBLUPsDF$DWpM + allBLUPsDF$PDW, allBLUPsDF[,-1])
# saveRDS(allBLUPsDF, file="allBLUPsDF_analyzeNH_20.rds")
# 
# # Best sporophytes
# bestSP_DWpM <- allBLUPsDF[progenySP,][order(allBLUPsDF$DWpM[progenySP], decreasing=T)[1:20],]
# bestSP_DWpM[,-(1:2)] <- bestSP_DWpM[,-(1:2)] %>% round(3)
# write.table(bestSP_DWpM, "BestSPbyDryWgtPerM_20.txt", quote=F, row.names=T, col.names=T)
# bestSP_idx <- allBLUPsDF[progenySP,][order(allBLUPsDF$index[progenySP], decreasing=T)[1:20],]
# bestSP_idx[,-(1:2)] <- bestSP_idx[,-(1:2)] %>% round(3)
# write.table(bestSP_idx, "BestSPbyDWpMandPDW_20.txt", quote=F, row.names=T, col.names=T)
# 
# # Best gametophytes
# bestGP_DWpM <- allBLUPsDF[releaseGP,][order(allBLUPsDF$DWpM[releaseGP], decreasing=T)[1:20],]
# bestGP_DWpM[,-(1:2)] <- bestGP_DWpM[,-(1:2)] %>% round(3)
# write.table(bestGP_DWpM[,-1], "BestGPbyDryWgtPerM_20.txt", quote=F, row.names=T, col.names=T)
# bestGP_idx <- allBLUPsDF[releaseGP,][order(allBLUPsDF$index[releaseGP], decreasing=T)[1:20],]
# bestGP_idx[,-(1:2)] <- bestGP_idx[,-(1:2)] %>% round(3)
# write.table(bestGP_idx[,-1], "BestGPbyDWpMandPDW_20.txt", quote=F, row.names=T, col.names=T)
# 
# # Ordered female and male gametophytes
# temp <- allBLUPsDF[allBLUPsDF$gpSex == "F",-(1:2)]
# temp <- temp[order(temp$DWpM, decreasing=T),] %>% round(3)
# write.table(temp, "FemaleGP_OrderedByDryWgtPerM_20.txt", quote=F, row.names=T, col.names=T)
# 
# temp <- allBLUPsDF[allBLUPsDF$gpSex == "M",-(1:2)]
# temp <- temp[order(temp$DWpM, decreasing=T),] %>% round(3)
# write.table(temp, "MaleGP_OrderedByDryWgtPerM_20.txt", quote=F, row.names=T, col.names=T)
# 
# 
# 


 














### Heritability conditional on development level at painting using the H matrix
# The checks don't have this development measure, so assign at random and
# do multiple times (nRepeats)

# nRepeats <- 2
# 
# nIsNAdev <- sum(is.na(dataNHpi$development))
# fracDevelIs1 <- sum(dataNHpi$development == 1, na.rm=T) / sum(!is.na(dataNHpi$development))
# develBLUPs <- develVu <- develVe <- h2covarPSi <- NULL
# 
# for (multImp in 1:nRepeats){
#   dataNHpi$develImp <- dataNHpi$development
#   dataNHpi$develImp[is.na(dataNHpi$development)] <- ifelse(runif(nIsNAdev) < fracDevelIs1, 1, 4)
#   msOutDevelh <- mixed.solve(y=dataNHpi$develImp, Z=msZ, K=hMat, X=msX, SE=T)
#   
#   develBLUPs <- cbind(develBLUPs, msOutDevelh$u)
#   develVu <- c(develVu, msOutDevelh$Vu)
#   develVe <- c(develVe, msOutDevelh$Ve)
#   
#   # Conditional heritabilities: Paint SP devel
#   # msXps <- model.matrix( ~ develImp + line%in%Year + block%in%Year+Year + popChk, data=dataNHpi) ### !!!???
#   
#   ###Year 2020 Only has 1 for develImp
#   msXps <- model.matrix( ~ develImp+line + block + popChk, data=dataNHpi) ### !!!???
#   
#   # All ES plots develImp is 1s for 2020 
#   
#   msOutWWPps <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msXps, SE=T)
#   
#   #Error in mixed.solve(y = log(dataNHpi$wetWgtPlot + 1), Z = msZ, K = hMat,  : X not full rank
#   
#   
#   msOutDWPMps <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msXps, SE=T)
#   msOutPDWps <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msXps, SE=T)
#   
#   h2covarPS <- c(heritability(msOutWWPps), heritability(msOutDWPMps), heritability(msOutPDWps))
#   names(h2covarPS) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")
#   h2covarPSi <- cbind(h2covarPSi, h2covarPS)
#   
# }

# msOutDevelh$u <- rowMeans(develBLUPs)
# msOutDevelh$Vu <- mean(develVu)
# msOutDevelh$Ve <- mean(develVe)
# h2covarPS <- rowMeans(h2covarPSi)



#  wetWgtPlot    dryWgtPerM    percDryWgt densityBlades  paintSPdevel 
#  0.369         0.357         0.014         0.219         0.512 

#
#

#  wetWgtPlot    dryWgtPerM    percDryWgt densityBlades  paintSPdevel 
#  0.470         0.468         0.008         0.134         0.307 




################################ ??????????????????????Need to fix this ???????!!!!!!!!!!!
# Function to calculate the mean, accounting for full rank incidence matrix
# WARNING: this function is hard-coded for the fixed effect incidence matrix
# msX <- model.matrix( ~ line + block + date + popChk, data=dataNHpi)
# So four lines, nine blocks, three dates

getMean <- function(beta){
  beta[1] + sum(beta[2:4])/4 + sum(beta[5:12])/9 + sum(beta[13:14])/3
}




intercepts <- c(getMean(msOutWWPh$beta), getMean(msOutDWPMh$beta), getMean(msOutPDWh$beta), getMean(msOutDBh$beta), getMean(msOutDevelh$beta))
names(intercepts) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt", "densityBlades", "paintSPdevel")
# wetWgtPlot and dryWgtPerM log transformed
intercepts[1:2] <- exp(intercepts[1:2]) - 1

# locMeansH <- msOutWWPh$beta[17:22]
# locMeansH <- locMeansH - min(locMeansH)
# locMeansSEh <- msOutWWPh$beta.SE[17:22]
# 

##################################





#####################

interceptsIM <- c(getMean(msOutBLh$beta), getMean(msOutBMWh$beta), getMean(msOutBW10h$beta), getMean(msOutBTh$beta), getMean(msOutSLh$beta), getMean(msOutSDh$beta))
names(interceptsIM) <- c("bladeLength", "bladeMaxWidth", "bladeWid10cm", "bladeThickness", "stipeLength", "stipeDiameter")
interceptsIM <- exp(interceptsIM) - 1

# Can't remember exactly why I wanted to do a multikernel.  I think having
# to do with the different levels on individuals and plots.  But I think
# the model did not converge
# emZlist <- list(emZsp, emZpl)
# emKlist <- list(hMat, diag(ncol(emZpl)))
# emOutBLh <- EMMREML::emmremlMultiKernel(y=log(dataNHim$bladeLength+1), X=msXim, Zlist=emZlist, Klist=emKlist)
# 



# Testing for significance of line and block
# ############################
# Here's how to fit the model
# `entryc` is the fixed effect.  It has a different value for each of the checks
# Say you had 3 checks, then you could have "A", "B", and "C" for the checks
# Plus one other unique value that is the _same_ for all the experimental entries
# Say you had 100 experimental entries, then, in addition to the A, B, and C for
# the checks, you could put in a value of 0 for each experimental entry
# So there would be 100 "0" in there.
# `blk` is just the levels of the blocks.  So if you have 20 blocks, it goes from 1 to 20.
# `entry` has a unique level for all different entries, be they checks or experimental
# So with the example of 3 checks and 100 experimental entries, it would go from 1 to 103
# `new` has two levels, one level for checks and one for experimental entries
# These can be 0 and 1, respectively
###############################




library(lme4)
exptlSP <- as.character(dataNHpi$popChk) == "ES"
dataNHpi$entry <- as.character(dataNHpi$popChk)
dataNHpi$entry[exptlSP] <- as.character(dataNHpi$plotNo[exptlSP])
dataNHpi$new <- as.factor(ifelse(exptlSP, 1, 0))


fitAug <- lmer(log(wetWgtPlot+1) ~ popChk + line*block + (1|entry:new), data=dataNHpi) #+ Year
print("Wet Weight Per Plot")
print(aov <- anova(fitAug))
# Significance tests when using the line * block interaction mean square as the error term
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff] / aov$'F value'[5], aov$npar[eff], aov$npar[5]), "\n")
# Significance tests when using the error mean square as the error term
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff], aov$npar[eff], 178), "\n")


print("Dry Weight Per Meter")
fitAug <- lmer(log(dryWgtPerM+1) ~ popChk + line*block + (1|entry:new), data=dataNHpi)
print(aov <- anova(fitAug))
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff] / aov$'F value'[5], aov$npar[eff], aov$npar[5]), "\n")
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff], aov$npar[eff], 178), "\n")
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff], aov$npar[eff], 24), "\n")


print("Percent Dry Weight")
fitAug <- lmer(percDryWgt ~ popChk + line*block  + (1|entry:new), data=dataNHpi)
print(aov <- anova(fitAug))
for (eff in 1:4) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff] / aov$'F value'[5], aov$npar[eff], aov$npar[5]), "\n")
for (eff in 1:4) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff], aov$npar[eff], 178), "\n")



# Measures of few kelp Protein Fat Fiber Ash
sum(!is.na(dataNHpi$protPerc))
pdf("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/Composition.pdf", height=5, width=10)
par(mfrow=c(1, 2))
boxplot(dataNHpi[, c("protPerc", "fiberPerc", "ashPerc")], names=c("Protein", "Fiber", "Ash"), ylim=c(5, 45), ylab="Percent", cex.lab=1.3, cex.axis=1.3)
boxplot(dataNHpi$fatPerc, ylim=c(0, 1.4), ylab="Percent", cex.lab=1.3)
mtext("Fat", side=1, line=1, cex=1.3)
dev.off()

# Filter to retain only the ten longest blades
# Some plots had many more blades measured, for some reason so down-sample at 
# random to 15, therefore do multiple times
nRep <- 5
h2onMI <- BLUPsOnMI <- list()
for (i in 1:nRep){
  imds <- dataNHim
  imds$plotNo <- as.character(imds$plotNo)
  # Down sample and pick ten with longest blades
  for (plot in unique(imds$plotNo)){
    # Down sample to fifteen blades
    whichPlot <- which(imds$plotNo == plot)
    nBlades <- length(whichPlot)
    if (nBlades > 15){
      imds <- imds[-sample(whichPlot, nBlades - 15),]
    }
    # Pick the ten biggest remaining
    whichPlot <- which(imds$plotNo == plot)
    nBlades <- length(whichPlot)
    if (nBlades > 10){
      blOrder <- order(imds$bladeLength[whichPlot])
      imds <- imds[-whichPlot[blOrder[1:(nBlades - 10)]],]
    }
  }
  
  # Complicated by the fact that there are no data for some im plots
  emZpl <- model.matrix( ~ -1 + plotNo, data=imds)
  imds$plotNo <- factor(imds$plotNo, levels=levels(dataNHpi$plotNo))
  msXim <- model.matrix( ~ -1 + plotNo, data=imds)
  msZim <- emZsp <- msXim %*% msZ
  msXim <- msXim %*% msXdb
  
  # bladeLength, bladeMaxWidth, bladeWid10cm, bladeThickness, stipeLength, stipeDiameter, stipeHollow, bullations, fouling
  msOutBLh <- mixed.solve(y=log(imds$bladeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)  ########## NOT FULL RANK !!!!!!!!!!!!!!!!!!!
  
  msOutBMWh <- mixed.solve(y=log(imds$bladeMaxWidth+1), Z=msZim, K=hMat, X=msXim, SE=T)
  msOutBW10h <- mixed.solve(y=log(imds$bladeWid10cm+1), Z=msZim, K=hMat, X=msXim, SE=T)
  
  cat("Replication", i, "\n")
  
  msOutBTh <- mixed.solve(y=log(imds$bladeThickness+1), Z=msZim, K=hMat, X=msXim, SE=T)
  msOutSLh <- mixed.solve(y=log(imds$stipeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
  msOutSDh <- mixed.solve(y=log(imds$stipeDiameter+1), Z=msZim, K=hMat, X=msXim, SE=T)
  
  h2covarHim2 <- c(heritability(msOutBLh), heritability(msOutBMWh), heritability(msOutBW10h), heritability(msOutBTh), heritability(msOutSLh), heritability(msOutSDh))
  names(h2covarHim2) <- c("bladeLength", "bladeMaxWidth", "bladeWid10cm", "bladeThickness", "stipeLength", "stipeDiameter")
  h2onMI <- c(h2onMI, list(h2covarHim2))
  
  allBLUPs2 <- cbind(msOutBLh$u, msOutBMWh$u, msOutBW10h$u, msOutBTh$u, msOutSLh$u, msOutSDh$u)
  colnames(allBLUPs2) <- c("BLen", "BMax", "B10", "BThk", "SLen", "SDia")
  BLUPsOnMI <- c(BLUPsOnMI, list(allBLUPs2))
}#END repeat with random sampling

meanh2onMI <- Reduce('+', h2onMI) / nRep
meanBLUPsOnMI <- Reduce('+', BLUPsOnMI) / nRep
allCor <- cor(meanBLUPsOnMI, allBLUPs[, -(1:5)])
print(diag(allCor))

pdf("~/Google Drive/Teaching/Seminars/2020/ARPA-E2020/Sources/CorrPlotLongTen.pdf")
corrplot::corrplot.mixed(cor(cbind(allBLUPs[,1:5], meanBLUPsOnMI)), diag="n", tl.cex=0.9, tl.col=1)
dev.off()

forHerit <- c(h2hMat, meanh2onMI)
rn <- which(names(forHerit)=="prePlantSPdevel")
if (length(rn) > 0) names(forHerit)[rn] <- "paintSPdevel"
rn <- which(names(forHerit)=="bladeWidth10cm")
if (length(rn) > 0) names(forHerit)[rn] <- "bladeWid10cm"

pdf("~/Google Drive/Teaching/Seminars/2020/ARPA-E2020/Sources/Heritabilities.pdf", height=6, width=10)
par(mar=par("mar")+c(4,0,0,0))
barplot(forHerit, ylab="Heritability", las=2, cex.axis=1.3, cex.lab=1.3, cex.names=1.3)
dev.off()

# Conditional heritabilities: Blade Density
msXdb <- model.matrix( ~ densityBlades + line + block + date + popChk, data=dataNHpi)

msOutWWPdb <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msXdb, SE=T)
msOutDWPMdb <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msXdb, SE=T)
msOutPDWdb <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msXdb, SE=T)

h2covarDB <- c(heritability(msOutWWPdb), heritability(msOutDWPMdb), heritability(msOutPDWdb))
names(h2covarDB) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")

forCondHerit <- c(forHerit[1], h2covarDB[1], h2covarPS[1], forHerit[2], h2covarDB[2], h2covarPS[2])
names(forCondHerit) <- c("", "wetWgtPlot", "", "", "dryWgtPerM", "")
pdf("~/Google Drive/Teaching/Seminars/2020/ARPA-E2020/Sources/CondHeritab.pdf", height=6, width=8)
# par(mar=par("mar")+c(2,0,0,0))
barplot(forCondHerit, ylab="Heritability", cex.axis=1.3, cex.lab=1.3, cex.names=1.3, col=c("gray", "dark green", "dark red"))
legend(2.57, 0.555, c("Conditional on Blade Density", "Conditional on Paint SP Devel."), pch=22, col=1, pt.bg=c("dark green", "dark red"), pt.cex=2, cex=1.5)
dev.off()


# # 2. Replace fixed effect of experimental versus check with the experimental location of origin
# # Both components will take out some genetic variation.  Should probably figure out how much each...
# 
# # Code to figure out incidence matrix of sampling location
# femaLocInc <- model.matrix(~ -1 + femaParLoc, data=dataNHpi)
# maleLocInc <- model.matrix(~ -1 + maleParLoc, data=dataNHpi)
# locInc <- (femaLocInc[,-1] + maleLocInc[,-1]) / 2
# # locInc sums to 1 for all experimental sporophytes but 0 for the checks
# # So perfectly colinear with the popChkES incidence column (17)
# # msXdb <- cbind(msX[, -17], locInc) 
# # 20190720 locInc causing a problem.  Not using it.