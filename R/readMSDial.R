  library(QFeatures)
  library(SummarizedExperiment)
  
  readMSDial <- function(file, version){
    table <- read.table(file, fill = TRUE, sep = "\t",
                        quote = "", header = FALSE)
    
    # Identify the starting row and column of the data
    startRow <- which(table[, 1] != "")[1]
    startCol <- which(table[1, ] != "")[1]
    ##TODO: version dependent error message if startRow or startCol are not as expected.
    
    # Split the table in parts
    colDataRaw <- table[1:startRow, startCol:ncol(table)]
    rowDataRaw <- table[startRow:nrow(table), 1:(startCol)]
    countsRaw <- table[startRow:nrow(table), startCol:ncol(table)]
    
    # Extract ids and counts data
    ids <- rowDataRaw[-1, 1]
    counts <- as.matrix(countsRaw[-1, -1])
    counts <- matrix(as.numeric(counts), nrow = nrow(counts), ncol = ncol(counts))
    colnames(counts) <- as.character(countsRaw[1, -1])
    rownames(counts) <- ids
    
    # Ensure row names of colData match counts column names
    colData <- DataFrame(t(colDataRaw[-nrow(colDataRaw), -1]))
    rownames(colData) <- as.character(colDataRaw[nrow(colDataRaw), -1])
    colnames(colData) <- as.character(colDataRaw[-nrow(colDataRaw), 1])

    # Ensure row names of rowData match counts row names
    rowData <- DataFrame(rowDataRaw[-1, ], row.names = ids)
    colnames(rowData) <- as.character(rowDataRaw[1,])

    
    
    # Create SummarizedExperiment object
   
    sumExp <- SummarizedExperiment(assays = list(counts = counts),
                                    rowData = rowData,
                                    colData = colData)
    ##TODO: Metadata with data source and version 
    
    # Create QFeatures object
    qf <- QFeatures()
    qf <- addAssay(qf, sumExp, name = "exampleAssay") 
    qf
    ##TODO: name
  }
  
 qf <- readMSDial("data/Metabolite_profile_showcase.txt")
 head(assay(qf))
colData(qf)