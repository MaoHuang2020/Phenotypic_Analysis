### Set working directory
# NOTE: Set working directory to that containing this code
# That directory is called "Code" here
# There should be a directory called "20200125" at the same level as "Code"
# which has the phenotypic data
# There should be a directory called "DArT" which has the genetic data
setwd("~/Documents/GitRepo/SugarKelpBreeding/TraitAnalyses200125/Code")
library(magrittr)
library(rrBLUP)

### Import data
# Plot level data
# spRows has 199 elements. dataNHpi has 218 rows.  The last 19 rows are checks
dataNHpi <- read.csv("../20200125/NH Phenotyping Data - Plot Info.csv", header=T, as.is=T)
dataNHpi <- dataNHpi[-(219:220),]
dataNHpi <- dataNHpi[,c(1, 3:6, 8:10, 14:15, 28:29, 42:44, 46:51)]
colnames(dataNHpi) <- c("plotNo", "femaPar", "femaParLoc", "malePar", "maleParLoc", "line", "block", "date", "wetWgtPlot", "lengthPlot", "numBlades5cm", "densityBlades", "percDryWgt", "wetWgtPerM", "dryWgtPerM", "development", "moistPerc", "protPerc", "fatPerc", "fiberPerc", "ashPerc")

# Data on individual blades
dataNHim <- read.csv("../20200125/NH Phenotyping Data - Individual Measurements.csv", header=T, as.is=T)
dataNHim <- dataNHim[, c(1, 3:12)]
colnames(dataNHim) <- c("plotNo", "subPlot", "bladeLength", "bladeMaxWidth", "bladeWid10cm", "bladeThickness", "stipeLength", "stipeDiameter", "stipeHollow", "bullations", "fouling")
dataNHim$plotNo <- gsub("S-", "", dataNHim$plotNo, fixed=TRUE)
dataNHim$plotNo[dataNHim$plotNo == ""] <- "C2-I"

# Marker data
source("convertDArTvcf.R")
fileName <- "../DArT/Report_DSacc18-3679_SNP_singlerow_2.csv"
fndrMrkData <- convertDArTvcf(fileName)
# Filtering of marker data
# Marker call rate > 0.8
nNAperMrk <- apply(fndrMrkData, 2, function(v) sum(is.na(v))) / nrow(fndrMrkData)
fndrMrkData <- fndrMrkData[, nNAperMrk < 0.2]
nNAperInd <- apply(fndrMrkData, 1, function(v) sum(is.na(v))) / ncol(fndrMrkData)
# Individual call rate > 0.8
fndrMrkData <- fndrMrkData[nNAperInd < 0.2,]
rownames(fndrMrkData) <- paste0(substring(rownames(fndrMrkData), first=1, last=2), "18", substring(rownames(fndrMrkData), first=3))
mrkRelMat <- A.mat(fndrMrkData, impute.method="EM", return.imputed=T, n.core=4)
fndrMrkDataImp <- mrkRelMat$imputed
mrkRelMat <- mrkRelMat$A

### Make pedigree relationship matrix
# Create the pedigree
source("makeBiphasicPed.R")
kelpNameColumns <- dataNHpi[, c("femaPar", "malePar")]
kelpNameColumns <- kelpNameColumns[kelpNameColumns[,1] != "",]
# rownames(mrkRelMat) gives the names of the founders to the pedigree
biphasicPedNH <- makeBiphasicPed(kelpNameColumns, rownames(mrkRelMat))

# Calculate the relationship matrix
source("calcCCmatrixBiphasic.R")
biphasicCCmat <- calcCCmatrixBiphasic(biphasicPedNH)
rownames(biphasicCCmat) <- colnames(biphasicCCmat) <- rownames(biphasicPedNH)
aMat <- 2 * biphasicCCmat
fndRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec == 0)))
gpRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) is.na(vec[2])))
spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))
nSp <- length(spRows)

# Combine the marker- with the pedigree- relationship matrix to make H matrix
hMat <- calcHmatrix(mrkRelMat, aMat, aMatFounders=rownames(mrkRelMat))
saveRDS(hMat, file="hMat_analyzeNH.rds")

###################################
## Data analysis really starts here

# Plot information analysis
# popChk: variable to indicate whether a plot is a check or an experimental sporophyte
dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "C", substr(dataNHpi$plotNo, 1, 2), "ES")
# withinLoc: variable is 1 if the cross was between gametophytes sampled from
# the same location, zero otherwise
dataNHpi$withinLoc <- c(ifelse(dataNHpi$femaParLoc[1:nSp] == dataNHpi$maleParLoc[1:nSp], 1, 0), rep(0, 19))

# Enforce data is numeric. If we don't have percent dry weight, set to missing
for (col in c("wetWgtPlot", "lengthPlot", "numBlades5cm", "densityBlades", "percDryWgt", "wetWgtPerM", "dryWgtPerM")) dataNHpi[,col] <- as.numeric(dataNHpi[,col])
dataNHpi[is.na(dataNHpi$percDryWgt), c("wetWgtPerM", "dryWgtPerM")] <- NA
keepRows <- !is.na(dataNHpi$percDryWgt)

# isSelf: variable is 1 if the cross was between gametophytes from the same 
# sporophyte, zero otherwise
isSelf <- integer(nrow(dataNHpi))
isSelf[diag(aMat[spRows, spRows]) > 1] <- 1
dataNHpi$isSelf <- isSelf

# Experimental design variables should be factors
for (col in c("plotNo", "femaPar", "femaParLoc", "malePar", "maleParLoc", "line", "block", "date", "popChk", "withinLoc", "isSelf")) dataNHpi[,col] <- as.factor(dataNHpi[,col])

msX <- model.matrix( ~ line + block + date + popChk, data=dataNHpi)
# Construction of Z assumes that all checks are at the bottom of the dataset
# Z has dimensions nrow(dataNHpi) x ncol(aMat)
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

### Heritability conditional on blade density using pedigree-based relationship
# Density Blades covariate
# The covariate relative to the focal trait could either remove sources of external random error or be a component that contributes genetic variation. The two are not mutually exclusive. To the extent that the former is prevalent the conditional heritability will be higher than the unconditional heritability. To the extent that the latter is prevalent, the reverse will hold. 

msXdb <- model.matrix( ~ densityBlades + line + block + date + popChk, data=dataNHpi)

msOutWWP <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=aMat, X=msXdb, SE=T)
msOutDWPM <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=aMat, X=msXdb, SE=T)
msOutPDW <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=aMat, X=msXdb, SE=T)

h2covar <- c(heritability(msOutWWP), heritability(msOutDWPM), heritability(msOutPDW))
names(h2covar) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")

locMeans <- msOutWWP$beta[17:22]
locMeans <- locMeans - min(locMeans)
locMeansSE <- msOutWWP$beta.SE[17:22]

### Heritability for traits using the H matrix
msOutWWPh <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msX, SE=T)
msOutDBh <- mixed.solve(y=dataNHpi$densityBlades, Z=msZ, K=hMat, X=msX, SE=T)
msOutDWPMh <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msX, SE=T)
msOutPDWh <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msX, SE=T)

### Heritability conditional on development level at painting using the H matrix
# The checks don't have this development measure, so assign at random and
# do multiple times (nRepeats)
nRepeats <- 10
nIsNAdev <- sum(is.na(dataNHpi$development))
fracDevelIs1 <- sum(dataNHpi$development == 1, na.rm=T) / sum(!is.na(dataNHpi$development))
develBLUPs <- develVu <- develVe <- h2covarPSi <- NULL
for (multImp in 1:nRepeats){
  dataNHpi$develImp <- dataNHpi$development
  dataNHpi$develImp[is.na(dataNHpi$development)] <- if_else(runif(nIsNAdev) < fracDevelIs1, 1, 4)
  msOutDevelh <- mixed.solve(y=dataNHpi$develImp, Z=msZ, K=hMat, X=msX, SE=T)
  develBLUPs <- cbind(develBLUPs, msOutDevelh$u)
  develVu <- c(develVu, msOutDevelh$Vu)
  develVe <- c(develVe, msOutDevelh$Ve)
  
  # Conditional heritabilities: Paint SP devel
  msXps <- model.matrix( ~ develImp + line + block + date + popChk, data=dataNHpi)
  
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

h2hMat <- c(heritability(msOutWWPh), heritability(msOutDWPMh), heritability(msOutPDWh), heritability(msOutDBh), heritability(msOutDevelh))
names(h2hMat) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt", "densityBlades", "paintSPdevel")

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

locMeansH <- msOutWWPh$beta[17:22]
locMeansH <- locMeansH - min(locMeansH)
locMeansSEh <- msOutWWPh$beta.SE[17:22]

##########################
# Inidividual measurements
# Somewhere here I would have to filter dataNHim to retain only the ten with
# longest blades
dataNHim$plotSub <- paste(dataNHim$plotNo, dataNHim$subPlot, sep="_")

for (col in c("bladeLength", "bladeMaxWidth", "bladeWid10cm", "bladeThickness", "stipeLength", "stipeDiameter")) dataNHim[,col] <- as.numeric(dataNHim[,col])

for (col in c("plotNo", "subPlot", "plotSub")) dataNHim[,col] <- as.factor(dataNHim[,col])

# Complicated by the fact that there are no data for some im plots
emZpl <- model.matrix( ~ -1 + plotNo, data=dataNHim)
levels(dataNHim$plotNo) <- levels(dataNHpi$plotNo)
msXim <- model.matrix( ~ -1 + plotNo, data=dataNHim)
msZim <- emZsp <- msXim %*% msZ
msXim <- msXim %*% msXdb

# bladeLength, bladeMaxWidth, bladeWid10cm, bladeThickness, stipeLength, stipeDiameter, stipeHollow, bullations, fouling
msOutBLh <- mixed.solve(y=log(dataNHim$bladeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBMWh <- mixed.solve(y=log(dataNHim$bladeMaxWidth+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBW10h <- mixed.solve(y=log(dataNHim$bladeWid10cm+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBTh <- mixed.solve(y=log(dataNHim$bladeThickness+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSLh <- mixed.solve(y=log(dataNHim$stipeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSDh <- mixed.solve(y=log(dataNHim$stipeDiameter+1), Z=msZim, K=hMat, X=msXim, SE=T)
h2hMatim <- c(heritability(msOutBLh), heritability(msOutBMWh), heritability(msOutBW10h), heritability(msOutBTh), heritability(msOutSLh), heritability(msOutSDh))
names(h2hMatim) <- c("bladeLength", "bladeMaxWidth", "bladeWid10cm", "bladeThickness", "stipeLength", "stipeDiameter")

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
names(h2hMat) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt", "densityBlades", "paintSPdevel")
names(h2hMatim) <- c("bladeLength", "bladeMaxWidth", "bladeWid10cm", "bladeThickness", "stipeLength", "stipeDiameter")

allBLUPs <- cbind(msOutWWPh$u, msOutDWPMh$u, msOutPDWh$u, msOutDBh$u, msOutDevelh$u, msOutBLh$u, msOutBMWh$u, msOutBW10h$u, msOutBTh$u, msOutSLh$u, msOutSDh$u)
colnames(allBLUPs) <- c("WWP", "DWpM", "PDW", "BDns", "SPdv", "BLen", "BMax", "B10", "BThk", "SLen", "SDia")

pdf("~/Google Drive/Teaching/Seminars/2020/ARPA-E2020/Sources/CorrPlot.pdf")
corrplot::corrplot.mixed(cor(allBLUPs), diag="n", tl.cex=0.9, tl.col=1)
dev.off()

sampledSP <- which(nchar(rownames(allBLUPs)) < 11)
releaseGP <- nchar(rownames(allBLUPs))
releaseGP <- which(10 < releaseGP & releaseGP < 15)
progenySP <- which(nchar(rownames(allBLUPs)) > 15)

fndNames <- rownames(allBLUPs)[sampledSP]
hasDesc <- sapply(fndNames, function(n) grep(n, rownames(allBLUPs)[progenySP]))
nDesc <- sapply(hasDesc, length)
colSPGP <- c(rep("black", length(sampledSP)), rep("dark green", length(releaseGP)), rep("dark red", length(progenySP)))
colSPGP[which(nDesc == 0)] <- "grey"
plot(allBLUPs[,2], pch=16, col=colSPGP, xlab="Individual number", ylab="Dry weight per plot")
rnPos <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]>0))]
rnNeg <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]<0))]
locPos <- table(substring(rnPos, 6, 7))
locNeg <- table(substring(rnNeg, 6, 7))

allBLUPsDF <- as.data.frame(allBLUPs)
allBLUPsDF <- cbind(pedigree=rownames(allBLUPs), allBLUPsDF)
rownames(allBLUPsDF)[progenySP] <- paste("SL18-UCONN-S", dataNHpi$plotNo[1:199], sep="") 

pdf("DWpMvsPDWsporophytes.pdf")
plot(allBLUPsDF$PDW[progenySP], allBLUPsDF$DWpM[progenySP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Sporophytes", cex.lab=1.3)
dev.off()
pdf("DWpMvsPDWgametophytes.pdf")
plot(allBLUPsDF$PDW[releaseGP], allBLUPsDF$DWpM[releaseGP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Gametophytes", cex.lab=1.3)
dev.off()

# Add gametophyte sex as a column to the dataframe
getGPsex <- function(gpName){
  sexNum <- strsplit(gpName, "-")[[1]][4]
  sex <- substring(sexNum, 1, 1)
  return(ifelse(sex %in% c("F", "M"), sex, NA))
}
gpSex <- sapply(rownames(allBLUPsDF), getGPsex)

# A selection index for sporophytes and gametophytes that weights DWpM twice as high as PDW
allBLUPsDF <- cbind(pedigree=allBLUPsDF$pedigree, gpSex=gpSex, index=2*allBLUPsDF$DWpM + allBLUPsDF$PDW, allBLUPsDF[,-1])
saveRDS(allBLUPsDF, file="allBLUPsDF_analyzeNH.rds")

# Best sporophytes
bestSP_DWpM <- allBLUPsDF[progenySP,][order(allBLUPsDF$DWpM[progenySP], decreasing=T)[1:20],]
bestSP_DWpM[,-(1:2)] <- bestSP_DWpM[,-(1:2)] %>% round(3)
write.table(bestSP_DWpM, "BestSPbyDryWgtPerM.txt", quote=F, row.names=T, col.names=T)
bestSP_idx <- allBLUPsDF[progenySP,][order(allBLUPsDF$index[progenySP], decreasing=T)[1:20],]
bestSP_idx[,-(1:2)] <- bestSP_idx[,-(1:2)] %>% round(3)
write.table(bestSP_idx, "BestSPbyDWpMandPDW.txt", quote=F, row.names=T, col.names=T)

# Best gametophytes
bestGP_DWpM <- allBLUPsDF[releaseGP,][order(allBLUPsDF$DWpM[releaseGP], decreasing=T)[1:20],]
bestGP_DWpM[,-(1:2)] <- bestGP_DWpM[,-(1:2)] %>% round(3)
write.table(bestGP_DWpM[,-1], "BestGPbyDryWgtPerM.txt", quote=F, row.names=T, col.names=T)
bestGP_idx <- allBLUPsDF[releaseGP,][order(allBLUPsDF$index[releaseGP], decreasing=T)[1:20],]
bestGP_idx[,-(1:2)] <- bestGP_idx[,-(1:2)] %>% round(3)
write.table(bestGP_idx[,-1], "BestGPbyDWpMandPDW.txt", quote=F, row.names=T, col.names=T)

# Ordered female and male gametophytes
temp <- allBLUPsDF[allBLUPsDF$gpSex == "F",-(1:2)]
temp <- temp[order(temp$DWpM, decreasing=T),] %>% round(3)
write.table(temp, "FemaleGP_OrderedByDryWgtPerM.txt", quote=F, row.names=T, col.names=T)

temp <- allBLUPsDF[allBLUPsDF$gpSex == "M",-(1:2)]
temp <- temp[order(temp$DWpM, decreasing=T),] %>% round(3)
write.table(temp, "MaleGP_OrderedByDryWgtPerM.txt", quote=F, row.names=T, col.names=T)

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
print("Wet Weight Per Plot")
fitAug <- lmer(log(wetWgtPlot+1) ~ popChk + line*block + date + (1|entry:new), data=dataNHpi)
print(aov <- anova(fitAug))
# Significance tests when using the line * block interaction mean square as the error term
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff] / aov$'F value'[5], aov$npar[eff], aov$npar[5]), "\n")
# Significance tests when using the error mean square as the error term
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff], aov$npar[eff], 178), "\n")

print("Dry Weight Per Meter")
fitAug <- lmer(log(dryWgtPerM+1) ~ popChk + line*block + date + (1|entry:new), data=dataNHpi)
print(aov <- anova(fitAug))
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff] / aov$'F value'[5], aov$npar[eff], aov$npar[5]), "\n")
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff], aov$npar[eff], 178), "\n")
for (eff in 1:5) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff], aov$npar[eff], 24), "\n")

print("Percent Dry Weight")
fitAug <- lmer(percDryWgt ~ popChk + line*block + date + (1|entry:new), data=dataNHpi)
print(aov <- anova(fitAug))
for (eff in 1:4) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff] / aov$'F value'[5], aov$npar[eff], aov$npar[5]), "\n")
for (eff in 1:4) cat(names(fitAug@frame)[eff+1], "\t", 1 - pf(aov$'F value'[eff], aov$npar[eff], 178), "\n")

# Measures of few kelp Protein Fat Fiber Ash
sum(!is.na(dataNHpi$protPerc))
pdf("~/Google Drive/Teaching/Seminars/2020/ARPA-E2020/Sources/Composition.pdf", height=5, width=10)
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
  msOutBLh <- mixed.solve(y=log(imds$bladeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
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

# Code to figure out incidence matrix of sampling location
femaLocInc <- model.matrix(~ -1 + femaParLoc, data=dataNHpi)
maleLocInc <- model.matrix(~ -1 + maleParLoc, data=dataNHpi)
locInc <- (femaLocInc[,-1] + maleLocInc[,-1]) / 2
# locInc sums to 1 for all experimental sporophytes but 0 for the checks
# So perfectly colinear with the popChkES incidence column (17)
# msXdb <- cbind(msX[, -17], locInc) 
# 20190720 locInc causing a problem.  Not using it.