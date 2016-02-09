sample_R <-
function (e, priR) {
  S <- t(e) %*% e + priR[2] * diag(ncol(e))
  v <- priR[1]
  sig <- solve(rwish(solve(S), v))
  cov2cor(sig)
}
