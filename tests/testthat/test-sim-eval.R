context("simulate & evaluate")

set.seed(1)

# create fake data
n <- 1E3
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

Y_sim = simulate(mock_bc, nsim = 1)[1, , ]

# Empirical probability of a community vector
p_Y = table(as.data.frame(Y)) / nrow(Y)
p_sim = table(as.data.frame(Y_sim)) / nrow(Y)


test_that("simulate.bayescomm works with no environment",{
  
  # How bad an approximation is p_sim for p_Y? Close to zero is good.
  kl_divergence = sum(p_Y * (log(p_Y) - log(p_sim)))
  
  expect_less_than(kl_divergence, .001)
})

test_that("likelihood evaluation works with no environment", {
  # Empirical log-likelihood from counting communities
  ll_empirical = sum(p_Y * log(p_Y))
  
  # log-likelihood from pmvnorm
  ll_bc = logLik(mock_bc, mock_bc$call$X[, -1, drop = FALSE], mock_bc$call$Y)
  
  expect_equal(ll_empirical, mean(ll_bc), tol = .01)
})
