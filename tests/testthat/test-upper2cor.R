context("upper2cor")

test_that("upper2cor works", {
  
  # Test with 2x2 matrix
  expect_equal(upper2cor(.1), matrix(c(1, .1, .1, 1), ncol = 2))
  
  R = cor(mtcars)
  
  # Test with a larger matrix
  expect_equal(
    upper2cor(R[upper.tri(R)]),
    R,
    check.attributes = FALSE
  )
})
