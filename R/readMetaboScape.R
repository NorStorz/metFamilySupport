library(readxl)
library(dplyr)
library(purrr)
library(QFeatures)
library(SummarizedExperiment)

#' Read Metaboscape Output File into a QFeatures Object
#'
#' This function reads a metabolite profile output file (.xlsx) from Metaboscape and 
#' converts it into a QFeatures object.
#'
#' @param file A character string specifying the path to the Metaboscape output file (Excel format).
#' @param version A character string specifying the version of Metaboscape used to generate the file.
#'   This parameter is currently not used.
#'
#' @return A QFeatures object containing:
#'   \itemize{
#'     \item An assay named "exampleAssay" with the metabolite counts.
#'     \item Row data (feature metadata) extracted from the input file.
#'     \item Column data (sample metadata) extracted from the sample names, including injection order and sample name.
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a Metaboscape output file named "data.xlsx":
#' qf <- readMetaboscape("data.xlsx")
#' 
#' # Examine the structure of the resulting QFeatures object
#' qf
#' 
#' # Access the assay data
#' assay(qf[["exampleAssay"]])
#' 
#' # Access the row data (feature metadata)
#' rowData(qf[["exampleAssay"]])
#' 
#' # Access the column data (sample metadata)
#' colData(qf)
#' }
#'
#' @importFrom readxl read_xlsx
#' @importFrom dplyr mutate across
#' @importFrom purrr map_dfc
#' @importFrom QFeatures QFeatures
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @details

#' @note 
#'
#' @seealso 
#' \code{\link[QFeatures]{QFeatures}} for more information on the QFeatures class.
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} for details on the underlying data structure.
#'
#' @references
#' 
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
  qf <- QFeatures(list(exampleAssay = sumExp), colData = colData(sumExp))
  qf
}

qf <- readMetaboscape("data/rye_test_data.xlsx")
a <- assay(qf)
