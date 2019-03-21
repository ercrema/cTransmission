#' @title Estimates the population frequencies of cultural variants.
#' @description The function estimates the population relative frequencies of cultural variants given observed raw frequencies and alpha priors.
#' @param x A vector of raw frequencies
#' @param nsim Number of repetitions. Default is 1. 
#' @param alpha Dirichlet \eqn{\alpha} parameter. Default is 1.
#' @return A matrix with relative frequencies of each variant (column) across \code{nsim} simulations (rows).
#' @export

generate_initialPop<-function(x, nsim=1, alpha=1)
{
  if (length(alpha)==1)
  {
    alpha = rep(alpha,length(x))
  }
  if(is.null(alpha)){
    stop("No alpha values defined!")
  } else{
    if (isTRUE(0 %in% alpha)){
      stop("All alpha values have to be larger than 0")
    }else{
      if(length(x)!=length(alpha)){
        stop("Length of alpha needs to be equal to the number of cultural variants + 1")
      }
    }
    exponents = x + alpha
  }
  return(rdirichlet(nsim,exponents))
}
