###GWAS on the GPSex
setwd("/Users/maohuang/Desktop/Kelp/SNP_calling/VarScan/SNPs_old_ref_genome/")
rm(list = ls())
inputData <- readRDS("inputData.rds")

mapData <- inputData$mapData[1:20,]
mrks <- inputData$mrks[,1:20]

  head(mrks)
  head(mapData)
phenos <- inputData$phenotypes
phenos$sexInd <- NA
phenos$sexInd[phenos$Sex == "MG"] <- 0
phenos$sexInd[phenos$Sex == "FG"] <- 1

mrkPs <- array()

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
GWASres$CHR <- "1"
GWASres$BP <- 1:nrow(GWASres)

colnames(GWASres)[3] <- "P"

pdf("ManhattanGamSex.pdf", h = 2.7, w = 5)
par(mar=c(3,3.5,1,1), mgp=c(1.8,0.5,0))
plot(GWASres$BP, -log(GWASres$P, 10), ylab = latex2exp::TeX("$-log_{10}(p)$"), xlab = "SNP index", cex = 0.5, pch = 19, col = "grey40", cex.lab = 0.8, cex.axis = 0.6, cex.main = 0.8)
dev.off()

