###### Founders SNPs
rm(list=ls())
library(data.table)
SNPs <- fread("fndrs_FilteredbiSNPs_numeric.raw", header = T, sep = "\t")
mapData <- read.table("fndrs_FilteredbiSNPs_numeric.map", header = F, sep = "\t")
#SNPs <- SNPs[1:(nrow(SNPs)-1) ,]
#phenos <- SNPs$IID
  dim(mapData)
SNPs.t<-t(SNPs[,-c(1:6)])
  dim(SNPs.t)
  dim(SNPs[,-c(1:6)])
  SNPs.t[1:5,1:6]
  SNPs[1:6,1:10]
colnames(SNPs.t)<-c(SNPs$IID)  
  SNPs.t[1:5,(ncol(SNPs.t)-5:ncol(SNPs.t))]
  mapData[1:5,]  
colnames(mapData)<-c("chrom","pos")  
  dim(mapData)
SNPs.t.m<-cbind(mapData,SNPs.t)
  dim(SNPs.t.m)  
  SNPs.t.m[1:5,1:6]
FounderSNP<-SNPs.t.m
save(FounderSNP,file="Founder_mainGenome.Rdata") 

#SNPs data
library(ggplot2)
Sites<-read.csv("SitesDP_hist.csv",sep=",",header=T)
head(Sites)
tail(Sites)
str(Sites)
#plot(Sites$Fraction_of_sites~Sites$BIN)
#levels(Sites$BIN)[2002]="2001"
colnames(Sites)<-c("Number_of_BIN","Fraction_of_sites")

dim(Sites)
Sites[2000:2002,]

max(Sites$Fraction_of_sites[-nrow(Sites)])
  #0.31
0.5*max(Sites$Fraction_of_sites[-nrow(Sites)])
  #0.156, choose 0.15 as the half of the max value
Sites[abs(Sites$Fraction_of_sites-0.156)<0.001,]
  #On both sides of the peack value, BIN=625 and797 Manually picked

ggplot(Sites,aes(x=Number_of_BIN,y=Fraction_of_sites))+
  geom_point()+
  geom_vline(aes(xintercept=791),color="red")

#geom_vline(aes(xintercept=633),color="red")+
Samples<-read.csv("SamplesStats_for_Filter1SNPs.csv",sep=",",header=T)

ggplot(Samples,aes(x=average.depth))+
  geom_histogram(binwidth=0.1)+
  geom_vline(aes(xintercept=0.5),color="red")
  
mean(Samples$average.depth)



########
setwd("/Users/maohuang/Desktop/Kelp/SNP_calling/VarScan/SNPs_old_ref_genome")
rm(list = ls())

library(data.table)
library(stringr)
library(ggplot2)

SNPs0 <- fread("SNPs_bcftoolsCalled_filter3_numeric.raw", header = T, sep = "\t")
#SNPs<-read.table("newFile_sub.txt", header = T, sep = "\t")
  dim(SNPs0)
  SNPs0[1:5,1:6]
  SNPs0[1:3,1:10]
  SNPs0[(nrow(SNPs0)-5:nrow(SNPs0)),1:10]
mapData0 <- read.table("SNPs_bcftoolsCalled_filter3_numeric.map", header = F, sep = "\t")
      #SNPs <- SNPs[1:(nrow(SNPs)-1) ,]
phenos0 <- SNPs0$IID
    head(phenos0)

# parse IID to extract phenotypes and covariates
phenos <- str_split_fixed(phenos0, "-", 5)
  head(phenos)
  dim(phenos)
phenos<-data.frame(phenos)
      #phenos <- phenos[,5:6]
      #phenos <- data.frame(cbind(str_split_fixed(phenos[,1], "-", 5), phenos[,2]))
colnames(phenos) <- c("Species", "Loc", "SPNum", "Sex", "GPNum")


SNPs <- as.matrix(SNPs0[,7:ncol(SNPs0)]); dim(SNPs)  ### The rownames of SNPs is the individual names
colnames(SNPs)<-as.vector(paste0(mapData0$V1,"_",mapData0$V2))
rownames(SNPs)<-phenos0
#paste0(phenos[,1],"_",phenos[,2],"_",phenos[,3],"_",phenos[,4],"_",phenos[,5])


##### Can do heterozygousity calculations here $$$###



SNPs[1:4,1:6]
SNPs[SNPs==1]=NA  ### Replace all Heter into NAs
    dim(SNPs)
#filter markers for missingness, the matrix  n xm
propMiss <- apply(SNPs, 2, function(x) sum(is.na(x))/length(x))
  str(propMiss)   #"000527F_317" 
  head(propMiss)
# count non-NAs in each row (marker) and non-NAs in each col(individual)
CallRateSNP<-apply(SNPs,2,function(x) sum(!is.na(x))/length(x))
  ##write.csv(CallRateSNP,"CallRateSNP.csv")

length(which(propMiss>0.3))
  missIndx <- which(propMiss > 0.3)  ## Further Filter markers with 50% of NAs, RMed  45232  markers! 996,516 -> 951284
  head(missIndx)                                  ## Filter markers with 30% of NAs, RMed 186780
  str(missIndx)
  
SNPs <- SNPs[,-missIndx]; dim(SNPs)
mapData <- mapData0[-missIndx,]  
    str(mapData)
    dim(SNPs)   #173 x 951284 if NAs 50%
                #173 x 764504 if NAs 30%
    
  #filter individuals for missingness with 75% NAs
  propMissInd <- apply(SNPs, 1, function(x) sum(is.na(x))/length(x))
    max(propMissInd)
    str(propMissInd)
  CallRateInd<-apply(SNPs,1,function(x) sum(!is.na(x))/length(x))  
  write.csv(CallRateInd,"CallRateInd.csv")
    max(CallRateInd)
    min(CallRateInd) 
    
missIndx <- which(propMissInd > 0.75)  ### RMed 90% NA individuals
    str(missIndx)

###!!!Caution, if this is going to be 0 elements, it will mess up the SNPs subsetting 
 #SNPs <- SNPs[-missIndx,]; dim(SNPs)
 #phenos <- phenos[-missIndx ,]


  #filter for allele frequency
 SNPs1<-SNPs
 
 SNPs<-SNPs1
 SNPs[SNPs == 2] <- 1
 af0 <- colMeans(SNPs, na.rm = T)     # mean per col= sum(1)/nrow(data_INdi)
  str(af0)
  head(af0)
 
af <- ifelse(af0 > 0.5, 1-af0, af0)   # if af>0.5, cal MAF 1-af
  
  write.csv(af0,"Allelic_Freq.csv")
  write.csv(af,"MAF.csv")
  MAF<-read.csv("MAF.csv",sep=",",header=T)
head(MAF)  
  colnames(MAF)<-c("Marker","MAF")
  ggplot(MAF,aes(x=MAF))+
    geom_histogram(binwidth=0.01)+
    geom_vline(xintercept=0.05,color="red")
    

Indx <- which(af < 0.05)
 ncol(SNPs)-length(afIndx)  
  str(afIndx) 
  head(afIndx)
  
SNPs <- SNPs[,!colnames(SNPs)%in%names(afIndx)]

#SNPs <- SNPs[,-afIndx]; dim(SNPs) # 536285 when from the Markers < 50% NAs  #470,960 SNPs when from the Markers <30% NAs
mapData1<-mapData 
mapData<-mapData1

mapData <- mapData[-afIndx,]
 
 saveRDS(list(phenotypes = phenos, mapData = mapData, mrks = SNPs), file="inputData.rds")

 
##################################
#### Assess the Heterozygousity on the SNPs(=SNPs0), which is 173 x996516 
 
HeterSNP0<-apply(SNPs, 2, function(x) length(which(x==1)))  ### Count the number of heterozygous
  head(HeterSNP0) 
  HeterSNP<-HeterSNP2
HeterSNP<-as.matrix(HeterSNP0)  
colnames(HeterSNP)<-"Hetercount"
  head(HeterSNP)
Contigs<- str_split_fixed(rownames(HeterSNP),"_",2)
Contigs<-data.frame(Contigs)
  head(Contigs)
  dim(Contigs)
  
HeterSNP3<-cbind(HeterSNP,Contigs)
Sum<-aggregate(HeterSNP3$Hetercount,list(Category=HeterSNP3$X1), FUN=sum)
  dim(Sum)
  head(Sum)
  
library(plyr)
Count<-count(HeterSNP3, "X1")
  dim(Count)  
  head(Count)
    identical(Count$X1,Sum$Category)
SumHet_NumSNP<-cbind(Sum,Count) 
  head(SumHet_NumSNP)
colnames(SumHet_NumSNP)<-c("Contig","Total_Number_of_Heter","Contig2","Number_of_SNPs")  
write.csv(SumHet_NumSNP[,-3],"Heterozygous_perContig.csv")   

HeterInd<-apply(SNPs, 1, function(x) length(which(x==1)))

Homo<-apply(SNPs,2,function(x)length(which(x==0 | x==2)))  

HomoInd<-apply(SNPs, 1, function(x) length(which(x==0 | x==2))) 

RatioInd<-HomoInd/HeterInd
  hist(RatioInd)
  head(HomoInd)
names(RatioInd)<-names(HomoInd)
  head(RatioInd)
write.csv(RatioInd,"Homo_to_Heter_RatioInd.csv")

    str(Homo)
    str(Heter)
  write.csv(Homo,"Homozygous.csv")
# Per individual
Ratio<-Homo/Heter 
  head(Ratio)
  write.csv(Ratio,"Homo_to_Heter_RatioSNPs.csv")

  Ratio<-read.csv("Homo_to_Heter.csv")
  ggplot(Rati,aes(x=HetNum))+
    geom_histogram(binwidth=1)
     
Heter<-read.csv("Heterozygous.csv",sep=",",header=T)

colnames(Heter)<-c("SNPID","HetNum")

  ggplot(Heter,aes(x=HetNum))+
    geom_histogram(binwidth=1)

##### Assess the CallRatemrk
CallRatemrk<-read.csv("CallRateSNP.csv",sep=",",header=T)
colnames(CallRatemrk)<-c("SNPID","CallRate")

  ggplot(CallRatemrk,aes(x=CallRate))+
    geom_histogram(binwidth=0.1)

##### Assess the CallRateInd
CallRateInd<-read.csv("CallRateInd.csv",sep=",",header=T)
colnames(CallRateInd)<-c("SampleID","CallRate")
  head(CallRateInd)

  ggplot(CallRateInd,aes(x=CallRate))+
    geom_histogram(binwidth=0.1)

min(CallRateInd$CallRate)
################################



rm(list = ls())
# saveRDS(list(phenotypes = phenos, mapData = mapData, mrks = SNPs), "analyses/inputData.rds")

# GWAS
inputData <- readRDS("inputData.rds")

mapData <- inputData$mapData
mrks <- inputData$mrks

# ### test a small set
# 
# mapData2<-mapData
# mapData<-mapData2[1:20,]
#   dim(mapData)
#   head(mapData)
#   
# mrks2<-mrks
# mrks<-mrks2[,1:20]
# ###


phenos <- inputData$phenotypes
phenos$sexInd <- NA
phenos$sexInd[phenos$Sex == "MG"] <- 0
phenos$sexInd[phenos$Sex == "FG"] <- 1

mrkPs <- array()

#cl <- makeCluster(8)
#registerDoParallel(cl)

for(m in 1:ncol(mrks)){
  tmpDF <- data.frame(sexInd = phenos$sexInd, Loc = phenos$Loc, mrk = mrks[,m])
  tmpDF <- tmpDF[!is.na(tmpDF$mrk) ,]
  fullmod <- glm(sexInd ~ mrk + Loc, data = tmpDF, family = "binomial")
  redmod <- glm(sexInd ~ Loc, data = tmpDF, family = "binomial")
  statistic <- 2 * (logLik(fullmod)[1] - logLik(redmod)[1])
  df1 <- redmod$df.residual
  df2 <- fullmod$df.residual
  df <- df1 - df2
  mrkPs[m] <- pchisq(statistic, df, lower.tail = FALSE)
}

GWASres <- mapData; GWASres$p <- mrkPs

inputData[["GWASres"]] <- GWASres

saveRDS(inputData, "inputData2.rds")

#Plot results. Since there are many many contigs (~3k) I plotted GWAS results using the SNP index rather than chromosome and position.

rm(list = ls())

inputData <- readRDS("inputData2.rds")
GWASres <- inputData$GWASres


# Order by scaffold
GWASres2<-GWASres[order(GWASres$V1,GWASres$V2),]

GWASres<-GWASres2
GWASres$CHR <- as.numeric(GWASres$V1)
colnames(GWASres)[3] <- "P"

GWASres$BP <- 1:nrow(GWASres)


pdf("GWAS_ManhattanGamSex.pdf", h = 2.7, w = 5)

tiff("GWAS_ManhattanGamSex.tiff",width=1200,height = 800,units = "px",pointsize=12)

par(mar=c(3,3.5,1,1), mgp=c(1.8,0.5,0))
plot(GWASres$BP, -log(GWASres$P, 10), ylab = latex2exp::TeX("$-log_{10}(p)$"), xlab = "SNP index", cex = 0.5, pch = 19, col = "grey40", cex.lab = 0.8, cex.axis = 0.6, cex.main = 0.8)

dev.off()

#
install.packages("qqman")
library(qqman)
GWASres$BP2<-GWASres$BP
GWASres$V1<-factor(GWASres$V1)

for (i in unique(levels(GWASres$V1))){
  GWASres[GWASres$V1==i,]$BP<-1:nrow(GWASres[GWASres$V1==i,])
}

write.csv(GWASres,"GWASres.csv")

GWASres$SNP<-paste0("SNP",1:nrow(GWASres))


tiff("GWAS_ManhattanGamSex_2.tiff",width=1200,height = 800,units = "px",pointsize=12)
par(mar=c(3,3.5,1,1), mgp=c(1.8,0.5,0))
manhattan(GWASres, chr="CHR", bp="BP", snp="SNP", p="P",ylim = c(0, 27), cex = 0.5,cex.lab=0.8,cex.axis=0.6,cex.main=0.8 )

dev.off()



####Count the number of SNPsby group
library(plyr)
NumSNPs<-count(GWASres, "CHR")
#OR
#library(qqman)
#NumSNPs<-as.data.frame(table(GWASres$CHR))

