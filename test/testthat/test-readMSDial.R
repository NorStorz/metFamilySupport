library(QFeatures)
source("R/readMsdial.R") 

qf <- readMSDial("data/Metabolite_profile_showcase.txt")
load("data/processMS1data.Rdata")
counts <- sapply(18:23,function(i) as.numeric(metaboliteProfile[,i]))


test_that("Sum of Intensities is correct", {
  count_cols <- grep("TRI|LVS", names(metaboliteProfile))
  counts <- sapply(count_cols,function(i) as.numeric(metaboliteProfile[,i]))
  
  sumQfeature <- sum(colSums(assay(qf)))
  sumMsdial <- sum(counts)
  expect_equal(sumQfeature, sumMsdial)

})

test_that("Rows and Columns are correct", {
  
  # check if number of rows is identical
  nrowQfeature <- nrow(assay(qf))
  nrowMsdial <- nrow(metaboliteProfile)
  expect_equal(nrowQfeature, nrowMsdial)
  
  # check if number of cols is identical
  ncolQfeature <- ncol(assay(qf))+ ncols(rowData(qf))
  ncolMsdial <- ncol(metaboliteProfile)
  expect_equal(as.integer(ncolQfeature), ncolMsdial)
  
  #check if colnames match
  colNamesQfeature <- c(colnames(rowData(qf))[[1]], rownames(colData(qf)))
  colNamesMsdial <- colnames(metaboliteProfile)
  
  l1 <- length(colNamesMsdial[!(colNamesMsdial %in% colNamesQfeature)])
  l2 <- length(colNamesQfeature[!(colNamesQfeature %in% colNamesMsdial)])
  
  expect_equal(l1, l2)
  expect_equal(l2, 0)
})


