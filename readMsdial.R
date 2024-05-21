library(QFeatures)
library(SummarizedExperiment)

readMSDial <- function(file, version){
  table <- read.table(file, fill = TRUE, sep = "\t",
                      quote = "", header = TRUE)
  
  # Identify the starting row and column of the data
  r <- which(table[, 1] != "")[1]
  c <- which(table[1, ] != "")[1]
  
  # Split the table in parts
  tr <- table[1:r, c:ncol(table)]
  bl <- table[r:nrow(table), 1:(c-1)]
  br <- table[r:nrow(table), c:ncol(table)]
  
  # Extract ids and matrix data
  ids <- bl[-1, 1]
  mat <- as.matrix(br[-1, -1])
  mat <- matrix(as.integer(mat), nrow = nrow(mat), ncol = ncol(mat))
  colnames(mat) <- as.character(br[1, -1])
  rownames(mat) <- ids
  
  # Ensure row names of rowdata_x match mat row names
  rowdata_x <- DataFrame(bl[-1, -1], row.names = ids)
  colnames(rowdata_x) <- as.character(bl[1, -1])
  
  # Ensure row names of coldata_x match mat column names
  coldata_x <- DataFrame(t(tr[-nrow(tr), -1]))
  rownames(coldata_x) <- as.character(tr[nrow(tr), -1])
  colnames(coldata_x) <- as.character(tr[-nrow(tr), 1])
  
  # Create SummarizedExperiment object
  test_se <- SummarizedExperiment(assays = list(counts = mat),
                                  rowData = rowdata_x,
                                  colData = coldata_x)
  
  # Create QFeatures object
  qf <- QFeatures()
  qf <- addAssay(qf, test_se, name = "exampleAssay")
  qf
}

qf <- readMSDial("Metabolite_profile_showcase.txt")
head(assay(qf))
