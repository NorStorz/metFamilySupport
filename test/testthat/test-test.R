test_that("add function works correctly", {
  add <- function(x, y) x + y
  expect_equal(add(1, 1), 2)
  expect_equal(add(2, 2), 4)
  expect_equal(add(-1, 1), 0)
})
