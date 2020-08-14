#' Calculate a coefficient of coancestry matrix
#'
#' \code{calcCCmatrixBiphasic} returns an additive relationship matrix from a
#'  pedigree specified in three columns. The first column has to be the row
#'  number and sire and dam columns refer directly to rows of the pedigree.
#'
#' \code{calcCCmatrixBiphasic} has functionality useful for species with a
#'  biphasic lifecycle. A haploid progeny can be specified by giving the
#'  row number of the diploid parent in the first column and setting the
#'  second column to NA. A diploid progeny should be specified by giving
#'  the row numbers of the two haploid parents.
#'
#'  The following rules are followed.
#'  1. Diploid from two diploids. The usual as given in quantitative genetics.
#'  2. Diploid from two haploids. The row is the mean of the haploid rows.
#'  The diagonal element is 0.5 * (1 + coancestry of the two haploids).
#'  3. Haploid from diploid. The row is the same as the row of the diploid.
#'  The diagonal element is 1.
#'  4. Haploid from haploid. The row is the same as the rwo of the haploid.
#'  The diagonal element is 1. (Equivalent to cloning the haploid).
#'
#' @param pedColumns A data.frame with three columns. The first column
#'  has to be the row number and sire and dam columns refer directly to rows
#'  of the pedigree. Parents of founders need to be set to 0 (ZERO). The row
#'  of a child has to be after (i.e. a higher row number) that of its parents.
#'  Indicate a haploid progeny by setting the second column to NA. If both
#'  columns have numbers in them, the progeny is assumed to be diploid.
#'  If an individual has one known and one unknown parent, set the unknown
#'  parent to 0.
#'
#' @return A matrix, \code{ccMat}, the coefficient of coancestry matrix
#'
calcCCmatrixBiphasic <- function(pedColumns){
    
    calcCCmatRow <- function(pedRec){ # Function to process one row of pedColumns
        prog <- pedRec[1]
        sire <- pedRec[2]
        dam <- pedRec[3]
        progM1 <- prog - 1
        if (sire){
            sireRow <- ccMat[sire, 1:progM1]
        } else{
            sireRow <- integer(progM1)
        }
        if (is.na(dam)){ # The progeny is a haploid
            ccMat[prog, 1:progM1] <<- sireRow
            ccMat[1:progM1, prog] <<- sireRow
            ccSelf <- 1
        } else{ # The progeny is a diploid
            if (dam){
                damRow <- ccMat[dam, 1:progM1]
            } else{
                damRow <- integer(progM1)
            }
            ccMat[prog, 1:progM1] <<- (sireRow + damRow) / 2
            ccMat[1:progM1, prog] <<- (sireRow + damRow) / 2
            if (sire > 0 & dam > 0){
                ccSelf <- (1 + ccMat[sire, dam]) / 2
            } else{
                ccSelf <- 0.5
            }
        }
        ccMat[prog, prog] <<- ccSelf
    }#END calcCCmatRow
    
    # calcCCmatRow <- function(pedRec){ # Function to process one row of pedColumns
    #     prog <- pedRec[1]
    #     sire <- pedRec[2]
    #     dam <- pedRec[3]
    #     progM1 <- prog - 1  # subtract 1, in the diagnol
    #     if (sire !=0){                          
    #         sireRow <- ccMat[sire, 1:progM1]  ### If sire !=0, working with Non-founder
    #     } else{
    #         sireRow <- integer(progM1)    ### Founder
    #     }
    #     if (is.na(dam)=TRUE){ # The progeny is a haploid, GP
    #         ccMat[prog, 1:progM1] <<- sireRow
    #         ccMat[1:progM1, prog] <<- sireRow
    #         ccSelf <- 1   # Diagonal
    #     } else{ # The progeny is a diploid, SP
    #         if (dam !=0){
    #             damRow <- ccMat[dam, 1:progM1]  # Non-founder, copy the row of the sire and the dam row
    #         } else{
    #             damRow <- integer(progM1)   # Founder , if Dam is founder it is 0
    #         } # CC is the average of those two
    #         ccMat[prog, 1:progM1] <<- (sireRow + damRow) / 2  
    #         ccMat[1:progM1, prog] <<- (sireRow + damRow) / 2
    #         if (sire > 0 & dam > 0){
    #             ccSelf <- (1 + ccMat[sire, dam]) / 2  # Figuring out itself
    #         } else{
    #             ccSelf <- 0.5
    #         }
    #     }
    #     ccMat[prog, prog] <<- ccSelf
    # }#END calcCCmatRow
    
    #pedColumns<-biphasicPedNH
    # calculate Coef Coan matrix here
    nInd <- nrow(pedColumns)
    # nInd<-5
    ccMat <- matrix(0, nInd, nInd)
    ccMat[1, 1] <- ifelse(is.na(pedColumns[1, 3]), 1, 0.5)
    apply(pedColumns[2:nInd,], 1, calcCCmatRow)
    rownames(ccMat) <- colnames(ccMat) <- pedColumns[,1]
    return(ccMat)
}

#' @param gMat A genomic relationship matrix.
#'
#' @param aMat A pedigree relationship matrix.
#'
#' @param aMatFounders The individuals to match in diagonal with gMat.
#'
#' @param relWgt Relative weights of the two matrices.
#'
#' @return A matrix, \code{ccMat}, the coefficient of coancestry matrix


#####################
#### 2016 matrix
calcHmatrix <- function(gMat, aMat, aMatFounders, wgtAmat=0.05){
    # Scale the gMat so its mean diagonal is equal to the mean diagonal of aMatFounders
    amfMean <- mean(diag(aMat[aMatFounders, aMatFounders]))
    gMatMean <- mean(diag(gMat))
    gMat <- gMat * amfMean / gMatMean
    # If there are individuals in gMat but not aMat, they will be made founders
    gNoA <- setdiff(rownames(gMat), rownames(aMat))
    if (length(gNoA) > 0){
        nGnoA <- length(gNoA)
        nAmat <- nrow(aMat)
        aMat <- cbind(rbind(diag(amfMean, nGnoA), matrix(0, nAmat, nGnoA)), rbind(matrix(0, nGnoA, nAmat), aMat))
        rownames(aMat)[1:nGnoA] <- colnames(aMat)[1:nGnoA] <- gNoA
    }
    # Make sure things are set up right to create the partition
    idxGinA <- sapply(rownames(gMat), function(n) which(rownames(aMat) == n))
    newOrder <- c(idxGinA, setdiff(1:nrow(aMat), idxGinA))
    aMat <- aMat[newOrder, newOrder]
    # Do the math
    idx11 <- 1:nrow(gMat)
    a11 <- aMat[idx11, idx11]
    a11inva12 <- solve(a11) %*% aMat[1:nrow(gMat), -(1:nrow(gMat))]
    gw <- (1 - wgtAmat)*gMat + wgtAmat*a11
    hMat <- cbind(rbind(gw, crossprod(a11inva12, gw)),
                rbind(gw, crossprod(a11inva12, gw - a11)) %*% a11inva12)  
    ### G_weighted, Cross Product, Transpose of the Cross product
    hMat[-idx11, -idx11] <- hMat[-idx11, -idx11] + aMat[-idx11, -idx11]
    rownames(hMat) <- colnames(hMat) <- rownames(aMat)
    return(hMat)
}
