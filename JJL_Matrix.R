
# Estimate heritabilities with more power (w/o response to E)
# Estimate line and block affect using checks (C-1 and C-2 etc.)
# Estimated breeding value 
# Pedigree tracking for measures of relatedness
# Use called SNPs for marker-based relationship matrix


rm(list=ls())

setwd("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/")
library(magrittr)
library(rrBLUP)

load("Dart2_mNA0.1_nNA0.5_RMhwe0.01_Manual_gl6m_0206.Rdata")
ls()
dim(gl6.m)
gl6.m[1:5,1:6]
Crossed<-as.vector(CrossedSP$Crossed)
str(Crossed)
gl6.m2<-gl6.m[rownames(gl6.m)%in%Crossed,]
dim(gl6.m2)
#write.csv(geno$taxa,"genotyped_Dart1.csv")

load("FarmCPU_GAPIT.Rdata")
geno2<-geno[,-1]
rownames(geno2)<-geno$taxa
geno2[1:5,1:6]

rownames(geno2) <- paste0(substring(rownames(geno2), first=1, last=2), "18", substring(rownames(geno2), first=3))
geno2[geno2==0]=-1
geno2[geno2==1]=0
geno2[geno2==2]=1

write.csv(rownames(geno2),"genotyped.csv")

##########This is JL marker scripting. 
# Output: -1, 0, 1 format matrix with individuals in rows and marker names in columns
# Marker data
# source("convertDArTvcf.R")
# fileName <- "../DArT/Report_DSacc18-3679_SNP_singlerow_2.csv"
# fndrMrkData <- convertDArTvcf(fileName)
# # Filtering of marker data
# # Marker call rate > 0.8
# nNAperMrk <- apply(fndrMrkData, 2, function(v) sum(is.na(v))) / nrow(fndrMrkData)
# fndrMrkData <- fndrMrkData[, nNAperMrk < 0.2]
# nNAperInd <- apply(fndrMrkData, 1, function(v) sum(is.na(v))) / ncol(fndrMrkData)
# # Individual call rate > 0.8
# fndrMrkData <- fndrMrkData[nNAperInd < 0.2,]
# rownames(fndrMrkData) <- paste0(substring(rownames(fndrMrkData), first=1, last=2), "18", substring(rownames(fndrMrkData), first=3))

fndrMrkData<-geno2
mrkRelMat <- A.mat(fndrMrkData, impute.method="EM", return.imputed=T)
fndrMrkDataImp <- mrkRelMat$imputed
mrkRelMat <- mrkRelMat$A


# Plot
dataNHpi <-read.csv("Plot_2019_2020_Compile_RM(dup)_CorrectName_FMG_SLSA18.csv",sep=",",header=TRUE)
  dim(dataNHpi)
  colnames(dataNHpi)

dataNHpi<-dataNHpi[!dataNHpi$crossID=="Buffer",]
#dataNHpi<-dataNHpi[!dataNHpi$crossID=="Check",]  ## RM checks????? # 569 plots

# Individual
dataNHim <- read.csv("Indi_2019_2020_Compile.csv",sep=",",header=TRUE)


dataNHim$plotNo <- gsub("S", "", dataNHim$plotNo, fixed=TRUE) ### The 2020 plotNo is "SXXX"; 2019 is just number
#dataNHim$plotNo[dataNHim$plotNo == ""] <- "C2-I"

tail(dataNHim)
#CrossedSP<-read.csv("Crossed_WildSP_2019_2020.csv",sep=",",header=TRUE)




### Make pedigree relationship matrix
# Create the pedigree
source("makeBiphasicPed.R")
kelpNameColumns <- dataNHpi[, c("femaPar", "malePar")]
kelpNameColumns <- kelpNameColumns[kelpNameColumns[,1] != "",]   #### ? Is this removing the Checks?


biphasicPedNH <- makeBiphasicPed(kelpNameColumns, rownames(mrkRelMat)) #### One has 18, the other not

        write.csv(biphasicPedNH,"biphasicPedNH.csv")
# Calculate the relationship matrix
source("calcCCmatrixBiphasic.R")
biphasicCCmat <- calcCCmatrixBiphasic(biphasicPedNH)

rownames(biphasicCCmat) <- colnames(biphasicCCmat) <- rownames(biphasicPedNH)
aMat <- 2 * biphasicCCmat
fndRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec == 0)))  ### Founders, which() only keeps the TRUE ones
gpRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) is.na(vec[2])))   ### GPs, the col 3 is NA
spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))    ### Progeny sps, col 2 and 3 are all values
nSp <- length(spRows)


# Combine the marker- with the pedigree- relationship matrix to make H matrix
hMat <- calcHmatrix(mrkRelMat, aMat, aMatFounders=rownames(mrkRelMat))
saveRDS(hMat, file="hMat_analyzeNH.rds")

length(fndRows) # 162
length(gpRows) #390
length(spRows) #569
