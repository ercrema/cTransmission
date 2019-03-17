#' @title Generate a random samples from the Dirichlet distribution
#' @param nsim Number of random vectors to generate
#' @param alpha Vector containing shape parameters.
#' @return A matrix containing randoms samples from the Dirichlet distribution, where each row is a one simulated set.
#' @export
rdirichlet<-function(nsim,alpha)
{
  n <- length(alpha)
  x <- matrix(rgamma(n * nsim, alpha), ncol = n , byrow = TRUE)
  res <- x %*% rep(1, n)
  x/as.vector(res)
}
