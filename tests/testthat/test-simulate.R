context("simulate")

test_that("simulate.bayescomm works",{
  set.seed(1)
  
  # create fake data
  n <- 1E4
  nsp <- 3
  
  W <- solve(rWishart(1, 10 + nsp, diag(nsp))[ , , 1])
  
  mu <- matrix(0, n, nsp)  # true mean
  p = pnorm(mu)
  
  e <- matrix(rnorm(n * nsp), n) %*% chol(W)  # true e
  z <- mu + e  # true z
  
  
  # Simulate with truncation ------------------------------------------------
  
  #  true presence/absence
  Y <- ifelse(z > 0, 1, 0)
  
  # Make a bayescomm object with two samples of R and B that match the true 
  # values
  mock_bc = list(
    trace = list(
      R = matrix(cov2cor(W)[upper.tri(W)], nrow = 2, ncol = 3, byrow = TRUE),
      B = replicate(nsp, list(matrix(0, nrow = 2, ncol = 2)))
    ),
    call = list(
      X = matrix(1, nrow = n, ncol = 2),
      Y = Y
    )
  )
  class(mock_bc) = "bayescomm"
  
  Y_sim = simulate(mock_bc, nsim = 1)
  
  
  p_Y = table(as.data.frame(Y)) / nrow(Y)
  p_sim = table(as.data.frame(Y_sim[1,,])) / nrow(Y)
  
  # How bad an approximation is p_sim for p_Y? Close to zero is good.
  kl_divergence = sum(p_Y * (log(p_Y) - log(p_sim)))
  
  expect_less_than(kl_divergence, .001)
})
