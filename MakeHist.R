# Make hist on the plots
# Record the ones has the sorus tissue

# Find the ones that are the same cross between both years based on the Crosses INFO?????

plotBLUPs<-read.csv("allBLUPsDF_BOTH_Plots.csv",sep=",",header=T)
  head(plotBLUPs)
library(ggplot2)  


histo<-ggplot(plotBLUPs,aes(x=index))+
  geom_histogram(binwidth=0.01)

segment_data = data.frame(
  x = c(plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"]),
  xend = c(plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"]), 
  y = rep(0,length(plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"])),
  yend = rep(5, length(plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"]))
)

 histo +
  geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="red")+
   labs(x="Index",y="Plot Count")
  
#### PDW  
 histo<-ggplot(plotBLUPs,aes(x=PDW))+
   geom_histogram(binwidth=0.01)
 
 segment_data = data.frame(
   x = c(plotBLUPs$PDW[plotBLUPs$Isolate_Sorus=="Yes"]),
   xend = c(plotBLUPs$PDW[plotBLUPs$Isolate_Sorus=="Yes"]), 
   y = rep(0,length(plotBLUPs$PDW[plotBLUPs$Isolate_Sorus=="Yes"])),
   yend = rep(10, length(plotBLUPs$PDW[plotBLUPs$Isolate_Sorus=="Yes"]))
 )
 
 histo +
   geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="red")+
   labs(x="PDW",y="Plot Count")

 #### DWpM  
 histo<-ggplot(plotBLUPs,aes(x=DWpM))+
   geom_histogram(binwidth=0.01)
 
 segment_data = data.frame(
   x = c(plotBLUPs$DWpM[plotBLUPs$Isolate_Sorus=="Yes"]),
   xend = c(plotBLUPs$DWpM[plotBLUPs$Isolate_Sorus=="Yes"]), 
   y = rep(0,length(plotBLUPs$DWpM[plotBLUPs$Isolate_Sorus=="Yes"])),
   yend = rep(5, length(plotBLUPs$DWpM[plotBLUPs$Isolate_Sorus=="Yes"]))
 )
 
 histo +
   geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="red")+
   labs(x="DWpM",y="Plot Count")
 
 
# geom_vline(xintercept=c(plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"]),linetype="dashed",color="red")+
# geom_abline(xintercept=plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"],col="red")
#hist(plotBLUPs$index)
#abline(v=c(plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"]),col="red")

 h2PlIn<-read.csv("Heritabilities_Plot_Indi.csv",sep=",",header=T)
 
 #wide to long format
 library(tidyr)
# Gather different traits,what is their value called, from col1:colN 
 h2PI<-gather(h2PlIn,Trait,Herit,wetWgtPlot:stipeDiameter,factor_key = T)
 head(h2PI)

levels(h2PI$Trait) <-c("WWP","DWpM","pDW","BL","BmWid","BTh","SL","SDia")
head(h2PI)

ggplot(h2PI,aes(fill=Data,y=Herit,x=Trait))+
      geom_bar(position=position_dodge(),stat="identity")
      

# 
h2PL<-read.csv("Heritabilities_Plot_Scenarios.csv",sep=",",header=T)
  head(h2PL) 

h2PL.long<-gather(h2PL,Trait,Herit,densityBlades:percDryWgt,factor_key=T)   
  head(h2PL.long)

##compared Amat vs Hmat (both using Model2, where densityBlades are conditioned)
h2AH<-h2PL.long[h2PL.long$Model==2,]
h2AH<-h2AH[!is.na(h2AH$Herit),]
h2AH

ggplot(h2AH,aes(fill=Matrix,y=Herit,x=Trait))+
  geom_bar(position=position_dodge(),stat="identity")+
    facet_grid(rows=vars(Data))
## compared model2 (condition on blade density) and model3 (condition on development stage)

h2M23<-h2PL.long[h2PL.long$Scenario=="c" |h2PL.long$Scenario=="d",]
h2M23<-h2M23[!h2M23$Trait=="densityBlades",]
  h2M23
h2M23$Scenario<-factor(h2M23$Scenario)
  str(h2M23)
levels(h2M23$Scenario)<-c("BladeDensity covariate","Development covariate")

ggplot(h2M23,aes(fill=Scenario,y=Herit,x=Trait))+
  geom_bar(position=position_dodge(),stat="identity")
