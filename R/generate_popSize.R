#' @title Estimate Population Sizes
#' @description Estimates population sizes from observed sample sizes given the sample fraction and interpolates across timesteps.
#' @param x a cFreqData class object
#' @param sMean mean of the ratio \eqn{s} between sample size and population size. 
#' @param sVariance variance of the ratio \eqn{s} between sample size and population size.
#' @details The function estimates population sizes for each timestep (from the first time-step to the start of the final time-step )
#' @return A vector containing the population size for each time-step.
#' @export


generate_popSize<-function(x, sMean, sVariance)
{
  if (!any(class(x)%in%c("cFreqData")))
  {
    stop("x must be a cFreqData class object")
  }
  
  Nh = rep(0,x$n.phases) #vector containing the size of the population for each time point ti
  s = rep(0,x$n.phases) #vector containing the ratio between sample size and popluation size for each time point ti
  
   
  for (i in 1:x$n.phases){
    s[i] = rnorm(1,sMean,sVariance) # different draw for each phase?
    while(s[i]<0){s[i]=rnorm(1,sMean,sVariance)}
    Nh[i] = round(x$n[i]/s[i])
  }
 
 #generating population sizes intermediate time points 
  N = rep(0, (x$tstamp[x$n.phases]-x$tstamp[1])+1)
  N[1] = Nh[1]
  for (i in 1:(x$n.phases-1)){
    time_diff = x$tstamp[i+1]-x$tstamp[i]+1
    #Linear interpolation between N[i] and N[i+1]
    h = (approx(c(1:2),c(Nh[i],Nh[i+1]),n = time_diff))
    N[(x$tstamp[i]-x$tstamp[1]+2):(x$tstamp[i+1]-x$tstamp[1]+1)] = round(h$y[2:time_diff])
  }
  return(N)
}
