#library(QFeatures)
#source("R/readMetaboScape.R") 


qf <- readMetaboscape("data/rye_test_data.xlsx")

test_that("Sum of Intensities is correct", {
  sumQf <- sum(colSums(assay(qf)))
  expect_equal(sumQf,307716661)
})

test_that("Number of Rows and Columns are correct", {
  
  # check if number of rows is identical
  nrowQf <- nrow(assay(qf))
  expect_equal(nrowQf, 805L)
  
  # check if number of cols is identical
  ncolQf <- ncol(assay(qf))+ ncols(rowData(qf))
  expect_equal(as.integer(ncolQf), 73)
  
})


