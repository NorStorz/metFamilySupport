qfmetab <- readMetaboscape("data/rye_test_data.xlsx")
qfmsdial <- readMSDial("data/Metabolite_profile_showcase.txt")


colnames(rowData(qfmsdial)[[1]])
colnames(rowData(qfmetab)[[1]])


####################################################################################
## aligned spectra
parsePeakAbundanceMatrix <- function(filePeakMatrix, 
                                     doPrecursorDeisotoping, 
                                     mzDeviationInPPM_precursorDeisotoping, 
                                     mzDeviationAbsolute_precursorDeisotoping, 
                                     maximumRtDifference, 
                                     progress=FALSE)
{
  ## read file
  if(!is.na(progress))  if(progress)  incProgress(amount = 0.1, detail = paste("Parsing MS1 file content...", sep = "")) else print(paste("Parsing MS1 file content...", sep = ""))
  
  dataFrameAll <- read.table(filePeakMatrix, header=FALSE, sep = "\t", as.is=TRUE, quote = "\"", check.names = FALSE, comment.char = "")
  dataFrameAll1 <- read.table(filePeakMatrix, header=TRUE, sep = "\t", as.is=TRUE, quote = "\"", check.names = FALSE, comment.char = "")
  
  oldFormat <- max(which(dataFrameAll[1:5, 1] == "")) == 3
  header_rowNumber <- ifelse(test = oldFormat, yes = 4, no = 5)
  dataFrameHeader <- dataFrameAll[1:header_rowNumber, ]# weg
  
  `%notin%` <- Negate(`%in%`)
  TR<- c("Reference RT","Reference m/z","Comment",
         "Manually modified for quantification",
         "Total score","RT similarity","Average","Stdev")
  IN2<-which(unname(dataFrameHeader[4,]) %notin% TR) #in rownames bla
  dataFrameHeader1 <-dataFrameHeader[IN2]
  
  dataFrame <- dataFrameAll[(header_rowNumber + 1):nrow(dataFrameAll), ]
  dataFrame1 <- dataFrame[IN2]
  
  colnames(dataFrame) <- dataFrameHeader[header_rowNumber, ]
  colnames(dataFrame1) <- dataFrameHeader1[header_rowNumber, ]
  
  numberOfPrecursors <- nrow(dataFrame1)
  numberOfPrecursorsPrior <- numberOfPrecursors
  
  columnIndexEndOfAnnotation <- max(match(x = "Class", 
                                          table = dataFrameHeader1[1, ]), 
                                    na.rm = TRUE) #match mit class aus coldata
  
  if(ncol(dataFrame1) > columnIndexEndOfAnnotation){
    dataColumnStartEndIndeces <- c(columnIndexEndOfAnnotation + 1, ncol(dataFrame1))
    numberOfDataColumns <- dataColumnStartEndIndeces[[2]] - dataColumnStartEndIndeces[[1]] + 1
    dataColumnNames <- colnames(dataFrame1)[dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
    
    sampleClass          <- dataFrameHeader1[1, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader1)]
    sampleType           <- dataFrameHeader1[2, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader1)]
    sampleInjectionOrder <- dataFrameHeader1[3, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader1)]
    batchID              <- NULL
    if(!oldFormat)
      batchID            <- dataFrameHeader1[4, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader1)]
  } else {
    dataColumnStartEndIndeces <- NULL
    numberOfDataColumns <- 0
    dataColumnNames <- NULL
    
    sampleClass          <- NULL
    sampleType           <- NULL
    sampleInjectionOrder <- NULL
    batchID              <- NULL
  }
  
  commaNumbers <- sum(grepl(x = dataFrame1$"Average Mz", pattern = "^(\\d+,\\d+$)|(^\\d+$)"))
  decimalSeparatorIsComma <- commaNumbers == nrow(dataFrame1)
  if(decimalSeparatorIsComma){
    if(!is.null(dataFrame1$"Average Rt(min)"))     dataFrame1$"Average Rt(min)"     <- gsub(x = gsub(x = dataFrame1$"Average Rt(min)", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Average Mz"))          dataFrame1$"Average Mz"          <- gsub(x = gsub(x = dataFrame1$"Average Mz", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Fill %"))              dataFrame1$"Fill %"              <- gsub(x = gsub(x = dataFrame1$"Fill %", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Dot product"))         dataFrame1$"Dot product"         <- gsub(x = gsub(x = dataFrame1$"Dot product", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Reverse dot product")) dataFrame1$"Reverse dot product" <- gsub(x = gsub(x = dataFrame1$"Reverse dot product", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    if(!is.null(dataFrame1$"Fragment presence %")) dataFrame1$"Fragment presence %" <- gsub(x = gsub(x = dataFrame1$"Fragment presence %", pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
    
    ## replace -1 by 0
    if(numberOfDataColumns > 0) {
      for(colIdx in dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]){
        dataFrame1[ , colIdx] <- gsub(x = gsub(x = dataFrame1[ , colIdx], pattern = "\\.", replacement = ""), pattern = ",", replacement = ".")
      }
    }
  }
  
  ###################
  ## column formats
  if(!is.null(dataFrame1$"Average Rt(min)"))     dataFrame1$"Average Rt(min)"     <- as.numeric(dataFrame1$"Average Rt(min)")
  if(!is.null(dataFrame1$"Average Mz"))          dataFrame1$"Average Mz"          <- as.numeric(dataFrame1$"Average Mz")
  if(!is.null(dataFrame1$"Fill %"))              dataFrame1$"Fill %"              <- as.numeric(dataFrame1$"Fill %")
  if(!is.null(dataFrame1$"MS/MS included"))      dataFrame1$"MS/MS included"      <- as.logical(dataFrame1$"MS/MS included")
  if(!is.null(dataFrame1$"Dot product"))         dataFrame1$"Dot product"         <- as.numeric(dataFrame1$"Dot product")
  if(!is.null(dataFrame1$"Reverse dot product")) dataFrame1$"Reverse dot product" <- as.numeric(dataFrame1$"Reverse dot product")
  if(!is.null(dataFrame1$"Fragment presence %")) dataFrame1$"Fragment presence %" <- as.numeric(dataFrame1$"Fragment presence %")
  
  #####################
  ## sorted by m/z (needed for deisotoping)
  if(!is.null(dataFrame1$"Average Mz"))
    dataFrame1 <- dataFrame1[order(dataFrame1$"Average Mz"), ]
  
  ## replace -1 by 0
  if(numberOfDataColumns > 0){
    for(colIdx in dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]){
      dataFrame1[ , colIdx] <- as.numeric(dataFrame1[ , colIdx])
      if(!is.na(sum(dataFrame1[,colIdx] == -1)))
        dataFrame1[(dataFrame1[,colIdx] == -1),colIdx] <- 0
    }
  }
  
  ## deisotoping
  numberOfRemovedIsotopePeaks <- 0
  if(doPrecursorDeisotoping & !is.null(dataFrame1$"Average Mz")){
    if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Precursor deisotoping...", sep = "")) else print(paste("Precursor deisotoping...", sep = ""))
    distance13Cminus12C <- 1.0033548378
    
    ## mark isotope precursors
    precursorsToRemove <- vector(mode = "logical", length = numberOfPrecursors)
    
    if(numberOfDataColumns > 0){
      intensities <- dataFrame1[ , dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
      medians <- apply(X = as.matrix(intensities), MARGIN = 1, FUN = median)
    }
    
    for(precursorIdx in seq_len(numberOfPrecursors)){
      if((precursorIdx %% (as.integer(numberOfPrecursors/10))) == 0)
        if(!is.na(progress))  if(progress)  incProgress(amount = 0.0, detail = paste("Precursor deisotoping ", precursorIdx, " / ", numberOfPrecursors, sep = "")) else print(paste("Precursor deisotoping ", precursorIdx, " / ", numberOfPrecursors, sep = ""))
      
      mzError <- dataFrame1$"Average Mz"[[precursorIdx]] * mzDeviationInPPM_precursorDeisotoping / 1000000
      mzError <- max(mzError, mzDeviationAbsolute_precursorDeisotoping)
      
      ## RT difference <= maximumRtDifference
      validPrecursorsInRt <- abs(dataFrame1$"Average Rt(min)"[[precursorIdx]] - dataFrame1$"Average Rt(min)"[-precursorIdx]) <= maximumRtDifference
      
      ## MZ difference around 1.0033548378 (first isotope) or 1.0033548378 * 2 (second isotope)
      validPrecursorsInMz1 <- abs((dataFrame1$"Average Mz"[[precursorIdx]] - distance13Cminus12C * 1) - dataFrame1$"Average Mz"[-precursorIdx]) <= mzError
      validPrecursorsInMz2 <- abs((dataFrame1$"Average Mz"[[precursorIdx]] - distance13Cminus12C * 2) - dataFrame1$"Average Mz"[-precursorIdx]) <= mzError
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
    dataFrame1 <- dataFrame1[!precursorsToRemove, ]
    
    numberOfRemovedIsotopePeaks <- sum(precursorsToRemove)
    numberOfPrecursors <- nrow(dataFrame1)
  }
  
  if(!is.na(progress))  if(progress)  incProgress(amount = 0, detail = paste("Boxing...", sep = "")) else print(paste("Boxing...", sep = ""))
  returnObj <- list()
  returnObj$dataFrame1 <- dataFrame1
  
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

