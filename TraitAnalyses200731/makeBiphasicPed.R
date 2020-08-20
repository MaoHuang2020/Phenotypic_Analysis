#' Make a three-column pedigree based on a data.frame formated as dataNH
#'
#' @param kelpNameColumns A two-column data.frame with female parents in
#' column 1 and male parents in column 2.  These columns have names in the
#' specific format of dataNH.
#' Name format: SA18-CB-10-F1
#' Two letters for the species, two digits for the sampling year
#' dash two letters for the sampling location dash
#' Sporophyte sample number dash
#' Single letter male or female, single digit gametophyt number
#'
#' @return A three-column pedigree matrix with row number and parent ID
#' listed as row numbers
#'
library(stringr)
makeBiphasicPed <- function(kelpNameColumns, founderFirst=NULL){
    founder <- as.character(unlist(kelpNameColumns))
    founder<-strsplit(founder,split="-",fixed=T)
    founder <- sapply(founder, function(vec) paste(vec[1:3], collapse="-"))
    founder <- unique(founder)
    founder <- c(founderFirst, setdiff(founder, founderFirst))
    nFounders <- length(founder)
        nFounders
    ####Getting the female & male list of names     
    parFemale <- as.character(unique(kelpNameColumns[,1]))
    parMale <- as.character(unique(kelpNameColumns[,2]))
    founderParFema <- sapply(strsplit(parFemale, split="-", fixed=T), function(vec) paste(vec[1:3], collapse="-"))
    founderParMale <- sapply(strsplit(parMale, split="-", fixed=T), function(vec) paste(vec[1:3], collapse="-"))
    
    fndRow <- 1:nFounders
    names(fndRow) <- founder
    
    #### Order the female and male in the proper list
    femaRow <- nFounders + 1:length(parFemale)
    names(femaRow) <- parFemale
    
    maleRow <- nFounders + length(parFemale) + 1:length(parMale)
    names(maleRow) <- parMale
    
    #write.csv(parFemale,"parFemale.csv")
    #write.csv(parMale,"parMale.csv")
    ####?
    fndPed <- cbind(fndRow, 0, 0)
    femaPed <- cbind(femaRow, fndRow[founderParFema], NA)
    malePed <- cbind(maleRow, fndRow[founderParMale], NA)
    sporProg <- 1:nrow(kelpNameColumns) + nFounders + length(parFemale) + length(parMale)
    
    #### This piece is wrong IF WITHOUT the as.character() part, which then mixed up the searching numbers
    names(sporProg)<-paste0(kelpNameColumns[,1], "x", kelpNameColumns[,2]) ### MH add !!
    
    progPed <- cbind(sporProg, femaRow[as.character(kelpNameColumns[,1])], maleRow[as.character(kelpNameColumns[,2])]) ### Here the rows are 1 extrac
    
 
    
    allPed <- rbind(fndPed, femaPed, malePed, progPed)
    rownames(allPed) <- c(founder, parFemale, parMale, paste(kelpNameColumns[,1], kelpNameColumns[,2], sep="x"))
    return(allPed)
}
