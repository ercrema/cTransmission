#' @title Computes the number of samples to be removed at each transition.
#' @param cFreqData class object
#' @param N vector of population sizes (user defined or created using \code{generate_N})
#' @param r population replacement rate
#' @return A data.frame containing the number of variants to be removed (u) and added (v) to the population at each time step.
#' @export

generate_removalReplacement <- function(x, N, r)
{
  u = v = rep(0, (x$tstamp[x$n.phases]-x$tstamp[1])+1)

  for (j in 1:(x$n.phases - 1)){
    n_step = x$tstamp[j+1] - x$tstamp[j] 
    if (j==(x$n.phases-1)){
      n_step = n_step 
    }
    for (i in 1:n_step) {
      if (N[(x$tstamp[j] - x$tstamp[1]) + i] <= N[(x$tstamp[j] - x$tstamp[1]) + i + 1]) {
        u[(x$tstamp[j] - x$tstamp[1]) + i] = floor(r * N[(x$tstamp[j] - x$tstamp[1]) + i])
        v[(x$tstamp[j] - x$tstamp[1]) + i] = N[(x$tstamp[j] - x$tstamp[1]) + i + 1] - (N[(x$tstamp[j] - x$tstamp[1]) + i] - u[(x$tstamp[j] - x$tstamp[1]) + i])
      } else{
        u[(x$tstamp[j] - x$tstamp[1]) + i] = floor(r * N[(x$tstamp[j] - x$tstamp[1]) + i + 1]) + (N[(x$tstamp[j] - x$tstamp[1]) + i] - N[(x$tstamp[j] - x$tstamp[1]) + i + 1])
        v[(x$tstamp[j] - x$tstamp[1]) + i] = N[(x$tstamp[j] - x$tstamp[1]) + i + 1] - (N[(x$tstamp[j] - x$tstamp[1]) + i] - u[(x$tstamp[j] - x$tstamp[1]) + i])
      }
    }
  }

  h = data.frame(u=u,v=v)
  return(h)
}



