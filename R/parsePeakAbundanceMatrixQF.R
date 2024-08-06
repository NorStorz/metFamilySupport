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
  
  dataFrameAll <- read.table(filePeakMatrix, header=FALSE, sep = "\t", as.is=TRUE, quote = "\"", check.names = FALSE, comment.char = "")
  
  
  #oldFormat <- max(which(dataFrameAll[1:5, 1] == "")) == 3
  header_rowNumber <- ncol(colData(qf))+1 
  cols_to_exclude <- c("Reference RT","Reference m/z","Comment",
                       "Manually modified for quantification",
                       "Total score","RT similarity","Average","Stdev") 
  
  
  rawDataFrameHeader <- dataFrameAll[1:header_rowNumber, ]

  
  cols_to_keep <- which(!(rawDataFrameHeader[header_rowNumber, ] %in% cols_to_exclude))
  # which(!(colnames(rowdata) %in% cols_to_exclude)) # konkatenieren mit Assay colnames ?
  dataFrameHeader <-rawDataFrameHeader[cols_to_keep]
  
  
  #dataFrame <- dataFrameAll[(header_rowNumber + 1):nrow(dataFrameAll), cols_to_keep] #coldata mit assay ?
  
  #combine rowData with assay 
  dataFrame <- cbind(rowData(qf)[[1]] ,assay(qf))

  dataFrame <- df[,cols_to_keep]
  
  nrows <- header_rowNumber
  ncols <- ncols(rowData(qf))
  empty <- as.data.frame(matrix("", nrow = nrows, ncol = ncols))
  coldata <- t(as.data.frame(colData(qf)))
  rownames(empty) <- rownames(coldata)
  colnames(empty) <- colnames(rowData(qf[[1]]))
  
  #dataFrameHeader <- cbind(empty, coldata)
  
 
  colnames(dataFrame) <- colnames(rowData(qf[[1]]))
  ncol(rowData(qf)[[1]]) # TODO: colnames wie viele wann benennen?
  numberOfPrecursors <- nrow(dataFrame)
  numberOfPrecursorsPrior <- numberOfPrecursors 
  
  columnIndexEndOfAnnotation <- max(match(x = "Class", 
                                          table = dataFrameHeader[1, ]), 
                                    na.rm = TRUE)
  
  if(ncol(dataFrame) > columnIndexEndOfAnnotation){
    dataColumnStartEndIndeces <- c(columnIndexEndOfAnnotation + 1, ncol(dataFrame))
    numberOfDataColumns <- dataColumnStartEndIndeces[[2]] - dataColumnStartEndIndeces[[1]] + 1
    dataColumnNames <- colnames(dataFrame)[dataColumnStartEndIndeces[[1]]:dataColumnStartEndIndeces[[2]]]
    
    sampleClass          <- dataFrameHeader[1, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader)]
    sampleType           <- dataFrameHeader[2, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader)]
    sampleInjectionOrder <- dataFrameHeader[3, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader)]
    batchID              <- NULL
    if(!oldFormat)
      batchID            <- dataFrameHeader[4, (columnIndexEndOfAnnotation + 1):ncol(dataFrameHeader)]
  } else {
    dataColumnStartEndIndeces <- NULL
    numberOfDataColumns <- 0
    dataColumnNames <- NULL
    
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
