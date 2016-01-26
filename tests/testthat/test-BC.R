context("BC results match true values")

# (Based on code from `example(BayesComm)`)

set.seed(1)

# create fake data
n <- 1000
nsp <- 15
k <- 3

X <- matrix(c(rep(1, n), rnorm(n * k)), n)  # covariate matrix

W <- matrix(rnorm(nsp * nsp), nsp)
W <- W %*% t(W) / 2  # true covariance matrix

# Using correlation instead of covariance to avoid misspecification
R <- cov2cor(W)

B <- matrix(rnorm(nsp * (k + 1), 0, 0.5), nsp)  # true covariates

mu <- apply(B, 1, function(b, x) x %*% b, X)  # true mean
p = pnorm(mu)

e <- matrix(rnorm(n * nsp), n) %*% chol(R)  # true e
z <- mu + e  # true z

Y <- ifelse(z > 0, 1, 0)  # true presence/absence

# Fit a glm to Y with an offset equal to probit(p).
# There should be no evidence that adding an intercept or a different slope
# improves the relationship between Y and p
p_glm <- glm(
  c(Y) ~ c(qnorm(p)) + offset(c(qnorm(p))), 
  family = binomial(link = "probit")
)
p_values <- summary(p_glm)$coefficients[ , 4]
expect_true(all(p_values > .01)) # Expect large p-values




# run BC (after removing intercept column from design matrix)
m1 <- BC(Y, X[, -1], model = "full", its = 1000)



test_that("true parameters are recovered", {
  R_diff <- colMeans(m1$trace$R) - R[upper.tri(R)]
  expect_true(
    mean(R_diff)^2 / var(R[upper.tri(R)]) < .01
  )
  
  
  B_diff <- c(t(sapply(m1$trace$B, colMeans))) - c(B)
  expect_true(
    mean(R_diff)^2 / var(R[upper.tri(R)]) < .01
  )
})


test_that("true probabilities are recovered", {
  p_array <- predict(m1, X[, -1])
  
  marginal_p <- apply(p_array, 2, rowMeans)
  
  # If the prediction code is correct, there shouldn't be any improvement after
  # accounting for the offset
  prediction_lm <- lm(
    c(qnorm(p)) ~ c(qnorm(marginal_p)) + offset(c(qnorm(marginal_p)))
  )
  
  # expect large p-values, fail to reject
  expect_true(all(coef(summary(prediction_lm))[ , 4] > .01))
})
