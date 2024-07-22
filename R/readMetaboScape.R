library(readxl)
library(dplyr)
library(purrr)
library(QFeatures)
library(SummarizedExperiment)

file <- "data/rye test data all features.xlsx"

readMetaboscape <- function(file, version){
  
  table <- read_xlsx(file)
  
  
  colnames <- colnames(table)
  colIdsSamples <- grepl("\\d+$", colnames)
  
  startOfSamples <- which(colIdsSamples)[1]
  
  
  
  
  ###COUNTS
  ids <- table[[1]]
  countsRaw <- table[,colIdsSamples]
  countsNumeric <- countsRaw %>%
  mutate(across(everything(), ~ as.numeric(.)))
  counts <- as.matrix(countsNumeric)
  rownames(counts) <- ids
  
  ###ROWDATA
  rowData <- DataFrame(table[,!colIdsSamples])
  
  ###COLDATA
  
  sampleNames <- colnames(table[,colIdsSamples])
  colDataRaw <- sapply(sampleNames, function(x) {
    pos <- max(gregexpr("[^0-9]", x)[[1]])
    c(substr(x, 1, (pos - 1)), substr(x, pos + 1, nchar(x)))
  })
  
  "Injection order" <- as.numeric(colDataRaw[2,])
  result_df <- DataFrame(result)
  colnames(result_df) <- c("First Part", "Second Part")

  ### MERGE TO QFEATURES
  sumExp <- SummarizedExperiment(assays = list(counts = counts),
                                 rowData = rowData,
                                 colData = colData)

  qf <- QFeatures()
  qf <- addAssay(qf, sumExp, name = "exampleAssay") 
  qf
}

qf <- readMetaboscape("data/rye test data all features.xlsx")
head(assay(qf))
colData(qf)


