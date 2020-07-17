setwd("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/")

### This file is modified manually to remove excessive columns
Yr2020<-read.csv("Plot_2020 Phenotyping Datasheet_07142020Download.csv",sep=",",header=TRUE)
head(Yr2020)
colnames(Yr2020)
colnames(Yr2020)<-c("plotNo", "crossID","femaPar", "femaParLoc", "malePar", "maleParLoc", "Crosses","PlantingDens","Development","line", "position","block", "PhotoScore","Year", "wetWgtPlot", "lengthPlot",  "wetWgtPerM", "percDryWgt","dryWgtPerM","densityBladesPerM")

Yr2019<-read.csv("Plot_2019_NH Phenotyping Data_0714Download.csv",sep=",",header=TRUE)
colnames(Yr2019)
head(Yr2019)

length(colnames(Yr2020)) #20
colnames(Yr2019)[1:20]<-c("plotNo", "crossID","femaPar", "femaParLoc", "malePar", "maleParLoc", "Crosses","Development","PlantingDens","line", "position","block", "PhotoScore","Year", "wetWgtPlot", "lengthPlot",  "percDryWgt","wetWgtPerM", "dryWgtPerM","densityBladesPerM")

Yr2019Sub<-Yr2019[,colnames(Yr2019)%in%colnames(Yr2020)]
  dim(Yr2019Sub)
  head(Yr2019Sub)
all.equal(colnames(Yr2019Sub),colnames(Yr2020))

Yr2019Sort<-Yr2019Sub[,colnames(Yr2020)]
  head(Yr2019Sort)
  head(Yr2019Sub)
  head(Yr2020)  
all.equal(colnames(Yr2019Sort),colnames(Yr2020))

Plot_Compile<-rbind(Yr2019Sort,Yr2020)  
  head(Plot_Compile)
  dim(Yr2019Sort)  
  dim(Yr2020)  
  dim(Plot_Compile)  

write.csv(Plot_Compile,"Plot_2019_2020_Compile.csv")



#### Individual Level    
Yr2020.Ind<-read.csv("Individual_2020 Phenotyping Datasheet_07142020Download.csv",sep=",",header=TRUE)
  head(Yr2020.Ind)
colkeep<-c("plotNo", "crossID","Year","sampleNo","bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter","SPwetWgt","SPpercDryWgt", "holdfast","stipeHollow", "bullations", "fouling")
  length(colkeep)
colnames(Yr2020.Ind)[1:length(colkeep)]<-colkeep
  head(Yr2020.Ind)

Yr2019.Ind<-read.csv("Individual_2019_NH Phenotyping Data_0714Download.csv",sep=",",header=TRUE)
  head(Yr2019.Ind)
colnames(Yr2019.Ind)<-c("plotNo", "crossID","Year","Subplot","bladeLength", "bladeMaxWidth", "bladeWidth10cm","bladeThickness", "stipeLength", "stipeDiameter", "stipeHollow", "bullations", "fouling","Sorus","Other") 


colBothYears<-c("plotNo", "crossID","Year","bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter","stipeHollow", "bullations", "fouling")

Yr2020.Ind_keep<-Yr2020.Ind[,colnames(Yr2020.Ind)%in%colBothYears]
  head(Yr2020.Ind_keep)

Yr2019.Ind_keep<-Yr2019.Ind[,colnames(Yr2019.Ind)%in%colBothYears]
  head(Yr2019.Ind_keep)
Yr2019.Ind_keep.Sort<-Yr2019.Ind_keep[,colnames(Yr2020.Ind_keep)]
  head(Yr2019.Ind_keep.Sort)
all.equal(colnames(Yr2019.Ind_keep.Sort),colnames(Yr2020.Ind_keep))


Indi_Compile<-rbind(Yr2019.Ind_keep.Sort,Yr2020.Ind_keep)
  dim(Indi_Compile)
  dim(Yr2020.Ind_keep)
  dim(Yr2019.Ind_keep.Sort)

write.csv(Indi_Compile,"Indi_Compile.csv")
