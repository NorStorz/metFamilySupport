#' Title
#'
#' @param qfeatures Object of type QFeatures 
#' @param doPrecursorDeisotoping 
#' @param mzDeviationInPPM_precursorDeisotoping 
#' @param mzDeviationAbsolute_precursorDeisotoping 
#' @param maximumRtDifference 
#' @param progress boolean whether or not to show a progress bar
#'
#' @return
#' @export
#'
#' @examples
parsePeakAbundanceMatrix <- function(qfeatures, 
                                     doPrecursorDeisotoping, 
                                     mzDeviationInPPM_precursorDeisotoping, 
                                     mzDeviationAbsolute_precursorDeisotoping, 
                                     maximumRtDifference, 
                                     progress=FALSE)
{
  ## read file

  if(!is.na(progress)) {
    if(progress) {
      incProgress(amount = 0.1, detail = paste("Parsing MS1 file content...", sep = ""))
    } else {
      print(paste("Parsing MS1 file content...", sep = ""))
    } 
  }
  
  cols_to_exclude <- c("Reference RT","Reference m/z","Comment",
                       "Manually modified for quantification",
                       "Total score","RT similarity","Average","Stdev") 
  
  cols_to_keep <- which(!colnames(rowData(qfmsdial))[[1]] %in% cols_to_exclude)

  dataFrame <- cbind(rowData(qfeatures)[[1]][,cols_to_keep] ,assay(qfeatures))

  ncol(rowData(qfeatures)[[1]]) 
  numberOfPrecursors <- nrow(dataFrame)
  numberOfPrecursorsPrior <- numberOfPrecursors 
  

  
  if(nrow(colData(qfeatures))>0){
    
    dataColumnStartEndIndeces <- 1
    numberOfDataColumns   <- nrow(colData(qfeatures))
    sampleClass           <- colData(qfeatures)$Class
    sampleType            <- colData(qfeatures)$Type
    sampleInjectionOrder  <- colData(qfmsdial)$"Injection order"
    batchID               <- NULL
    if(! is.null(colData(qfeatures)$BatchID))
      batchID            <- colData(qfeatures)$BatchID
    
    }   else {
    dataColumnStartEndIndeces <- NULL
    numberOfDataColumns <- 0
    sampleClass          <- NULL
    sampleType           <- NULL
    sampleInjectionOrder <- NULL
    batchID              <- NULL
  }
  
  commaNumbers <- sum(grepl(x = dataFrame$"Average Mz", pattern = "^(\\d+,\\d+$)|(^\\d+$)"))
  decimalSeparatorIsComma <- commaNumbers == nrow(dataFrame)
  if(decimalSeparatorIsComma){
    if(!is.null(dataFrame$"Average Rt(min)"))     dataFrame$"Average Rt(min)"     <- gsub(x = gsub(x = dataFrame$"Average Rt(min)", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame$"Average Mz"))          dataFrame$"Average Mz"          <- gsub(x = gsub(x = dataFrame$"Average Mz", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame$"Fill %"))              dataFrame$"Fill %"              <- gsub(x = gsub(x = dataFrame$"Fill %", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame$"Dot product"))         dataFrame$"Dot product"         <- gsub(x = gsub(x = dataFrame$"Dot product", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame$"Reverse dot product")) dataFrame$"Reverse dot product" <- gsub(x = gsub(x = dataFrame$"Reverse dot product", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame$"Fragment presence %")) dataFrame$"Fragment presence %" <- gsub(x = gsub(x = dataFrame$"Fragment presence %", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    
    ## replace -1 by 0
    if(numberOfDataColumns > 0) {
      for(colIdx in dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]){
        dataFrame[ , colIdx] <- gsub(x = gsub(x = dataFrame[ , colIdx], pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
      }
    }
  }
  
  ###################
  ## column formats
  if(!is.null(dataFrame$"Average Rt(min)"))     dataFrame$"Average Rt(min)"     <- as.numeric(dataFrame$"Average Rt(min)")
  if(!is.null(dataFrame$"Average Mz"))          dataFrame$"Average Mz"          <- as.numeric(dataFrame$"Average Mz")
  if(!is.null(dataFrame$"Fill %"))              dataFrame$"Fill %"              <- as.numeric(dataFrame$"Fill %")
  if(!is.null(dataFrame$"MS/MS included"))      dataFrame$"MS/MS included"      <- as.logical(dataFrame$"MS/MS included")
  if(!is.null(dataFrame$"Dot product"))         dataFrame$"Dot product"         <- as.numeric(dataFrame$"Dot product")
  if(!is.null(dataFrame$"Reverse dot product")) dataFrame$"Reverse dot product" <- as.numeric(dataFrame$"Reverse dot product")
  if(!is.null(dataFrame$"Fragment presence %")) dataFrame$"Fragment presence %" <- as.numeric(dataFrame$"Fragment presence %")
  
  #####################
  ## sorted by m/z (needed for deisotoping)
  if(!is.null(dataFrame$"Average Mz"))
    dataFrame <- dataFrame[order(dataFrame$"Average Mz"), ]
  
  ## replace -1 by 0
  if(numberOfDataColumns > 0){
    for(colIdx in dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]){
      dataFrame[ , colIdx] <- as.numeric(dataFrame[ , colIdx])
      if(!is.na(sum(dataFrame[,colIdx] == -1)))
        dataFrame[(dataFrame[,colIdx] == -1),colIdx] <- 0
    }
  }
  
  ## deisotoping
  numberOfRemovedIsotopePeaks <- 0
  if(doPrecursorDeisotoping & !is.null(dataFrame$"Average Mz")){
    if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Precursor deisotoping...", sep = "")) else print(paste("Precursor deisotoping...", sep = ""))
    distance13Cminus12C <- 1.0033548378
    
    ## mark isotope precursors
    precursorsToRemove <- vector(mode = "logical", length = numberOfPrecursors)
    
    if(numberOfDataColumns > 0){
      intensities <- dataFrame[ , dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
      medians <- apply(X = as.matrix(intensities), MARGIN = 1, FUN = median)
    }
    
    for(precursorIdx in seq_len(numberOfPrecursors)){
      if((precursorIdx %% (as.integer(numberOfPrecursors/10))) == 0)
        if(!is.na(progress))  if(progress)  incProgress(amount = 0.0, detail = paste("Precursor deisotoping ", precursorIdx, " / ", numberOfPrecursors, sep = "")) else print(paste("Precursor deisotoping ", precursorIdx, " / ", numberOfPrecursors, sep = ""))
      
      mzError <- dataFrame$"Average Mz"[[precursorIdx]] * mzDeviationInPPM_precursorDeisotoping / 1000000
      mzError <- max(mzError, mzDeviationAbsolute_precursorDeisotoping)
      
      ## RT difference <= maximumRtDifference
      validPrecursorsInRt <- abs(dataFrame$"Average Rt(min)"[[precursorIdx]] - dataFrame$"Average Rt(min)"[-precursorIdx]) <= maximumRtDifference
      
      ## MZ difference around 1.0033548378 (first isotope) or 1.0033548378 * 2 (second isotope)
      validPrecursorsInMz1 <- abs((dataFrame$"Average Mz"[[precursorIdx]] - distance13Cminus12C * 1) - dataFrame$"Average Mz"[-precursorIdx]) <= mzError
      validPrecursorsInMz2 <- abs((dataFrame$"Average Mz"[[precursorIdx]] - distance13Cminus12C * 2) - dataFrame$"Average Mz"[-precursorIdx]) <= mzError
      validPrecursorsInMz <- validPrecursorsInMz1 | validPrecursorsInMz2
      
      ## intensity gets smaller in the isotope spectrum
      if(numberOfDataColumns > 0){
        validPrecursorsInIntensity <- (medians[-precursorIdx] - medians[[precursorIdx]]) > 0
      } else {
        validPrecursorsInIntensity <- TRUE
      }
      
      if(any(validPrecursorsInRt & validPrecursorsInMz & validPrecursorsInIntensity))
        precursorsToRemove[[precursorIdx]] <- TRUE
    }
    
    ## remove isotopes
    dataFrame <- dataFrame[!precursorsToRemove, ]
    
    numberOfRemovedIsotopePeaks <- sum(precursorsToRemove)
    numberOfPrecursors <- nrow(dataFrame)
  }
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Boxing...", sep = "")) else print(paste("Boxing...", sep = ""))
  returnObj <- list()
  returnObj$dataFrame <- dataFrame
  
  ## meta
  returnObj$oldFormat <- oldFormat
  returnObj$numberOfPrecursors <- numberOfPrecursors
  returnObj$dataColumnStartEndIndeces <- dataColumnStartEndIndeces
  returnObj$numberOfDataColumns <- numberOfDataColumns
  
  ## group anno
  returnObj$sampleClass          <- sampleClass
  returnObj$sampleType           <- sampleType
  returnObj$sampleInjectionOrder <- sampleInjectionOrder
  returnObj$batchID              <- batchID
  
  ## misc
  returnObj$numberOfPrecursorsPrior <- numberOfPrecursorsPrior
  returnObj$numberOfRemovedIsotopePeaks <- numberOfRemovedIsotopePeaks
  
  return (returnObj)
}
