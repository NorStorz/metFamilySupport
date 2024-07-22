library(readxl)
library(dplyr)
library(purrr)
library(QFeatures)
library(SummarizedExperiment)


readMetaboscape <- function(file, version){
  
  table <- read_xlsx(file)
  
  
  colnames <- colnames(table)
  colIdsSamples <- grepl("\\d+$", colnames)
  
  startOfSamples <- which(colIdsSamples)[1]
  
  
  
  
  # Extract ids and counts data
  ids <- table[[1]]
  countsRaw <- table[,colIdsSamples]
  countsNumeric <- countsRaw %>%
  mutate(across(everything(), ~ as.numeric(.)))
  counts <- as.matrix(countsNumeric)
  rownames(counts) <- ids
  
  # Extract rowData
  rowData <- DataFrame(table[,!colIdsSamples])
  
  # Extract colData from sample Names
  sampleNames <- colnames(table[,colIdsSamples])
  colDataRaw <- sapply(sampleNames, function(x) {
    # find position of the first character before the Run number
    pos <- max(gregexpr("[^0-9]", x)[[1]])
      c(substr(x, 1, (pos - 1)), substr(x, pos + 1, nchar(x)))
  })
  colData <- DataFrame("Injection order" = colDataRaw[2,],
                       "Sample name" = colDataRaw[1,])
                    
  # Create SummarizedExperiment object
  sumExp <- SummarizedExperiment(assays = list(counts = counts),
                                 rowData = rowData,
                                 colData = colData)

  # Create QFeatures object
  qf <- QFeatures()
  qf <- addAssay(qf, sumExp, name = "exampleAssay") 
  qf
}

#qf <- readMetaboscape("data/rye test data all features.xlsx")