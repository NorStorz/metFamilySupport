library(QFeatures)
source("R/readMsdial.R") 

qf <- readMSDial("data/Metabolite_profile_showcase.txt")
load("data/processMS1data.Rdata")
counts <- sapply(18:23,function(i) as.numeric(metaboliteProfile[,i]))


test_that("Sum of Intensities is correct", {
  count_cols <- grep("TRI|LVS", names(metaboliteProfile))
  counts <- sapply(count_cols,function(i) as.numeric(metaboliteProfile[,i]))
  
  sum1 <- sum(colSums(assay(qf)))
  sum2 <- sum(counts)
  expect_equal(sum1, sum2)
  # Optional: Check for error handling
})

test_that("Rows and Columns are correct", {
  
  # check if number of rows is identical
  nrow1 <- nrow(assay(qf))
  nrow2 <- nrow(metaboliteProfile)
  expect_equal(nrow1, nrow2)
  
  # check if number of cols is identical
  ncol1 <- ncol(assay(qf))+ ncols(rowData(qf))
  ncol2 <- ncol(metaboliteProfile)
  expect_equal(as.integer(ncol1), ncol2)
  
  #check if colnames match
  colNames1 <- c(colnames(rowData(qf))[[1]], rownames(colData(qf)))
  colNames2 <- colnames(metaboliteProfile)
  
  l1 <- length(colNames2[!(colNames2 %in% colNames1)])
  l2 <- length(colNames1[!(colNames1 %in% colNames2)])
  
  expect_equal(l1, l2)
  expect_equal(l2, 0)
})
